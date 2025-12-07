import math

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    Var,
    Constraint,
    Expression,
    Objective,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.models.unit_models import (
    Product,
    Feed,
    Mixer,
    StateJunction,
    Separator,
    HeatExchanger,
    SplittingType,
    HeatExchangerFlowPattern,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.pressure_changer import Pump
from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap.costing import WaterTAPCosting

from srp.components.translator_sw_to_water import TranslatorSWtoWater
from srp.utils.utils import touch_flow_and_conc

__author__ = "Kurban Sitterley"

__all__ = [
    "build_bc",
    "set_bc_operating_conditions",
    "set_bc_scaling",
    "init_bc",
    "add_bc_costing",
    "scale_bc_costs",
    "run_sequential_decomposition",
    "solve_bc",
    "display_bc_flow_table",
    "report_pump",
    "report_bc",
    "print_bc_stream_flows",
]

solver = get_solver()

_log = idaeslog.getLogger("SRP")


def solve_bc(blk, tee=False):
    solver = get_solver()
    try:
        results = solver.solve(blk, tee=tee)
        print(f"termination MVC {results.solver.termination_condition}")
        assert_optimal_termination(results)
    except:
        results = solver.solve(blk, tee=tee)
        print(f"termination MVC {results.solver.termination_condition}")
        assert_optimal_termination(results)

    return results


def build_system(recovery=0.5, Qin=350, Cin=11408, feed_temp=27, **kwargs):

    Qin = Qin * pyunits.gallons / pyunits.minute
    Cin = Cin * pyunits.mg / pyunits.liter

    m = ConcreteModel()
    m.recovery_mass = recovery
    m.recovery_vol = recovery
    m.Qin = Qin
    m.Cin = Cin
    m.feed_temp = feed_temp

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.feed)
    m.fs.product = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.product)
    m.fs.disposal = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.disposal)

    m.fs.bc = bc = FlowsheetBlock(dynamic=False)

    build_bc(m, bc, **kwargs)

    m.fs.feed_to_mvc = Arc(source=m.fs.feed.outlet, destination=bc.feed.inlet)

    m.fs.mvc_to_product = Arc(source=bc.product.outlet, destination=m.fs.product.inlet)

    m.fs.mvc_to_disposal = Arc(
        source=bc.disposal.outlet, destination=m.fs.disposal.inlet
    )

    add_bc_costing(bc)

    m.fs.costing.base_currency = pyunits.USD_2023
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.product.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
    )

    m.fs.costing.heat_exchanger.material_factor_cost.fix(5)
    m.fs.costing.evaporator.material_factor_cost.fix(5)
    m.fs.costing.compressor.unit_cost.fix(1 * 7364)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_bc(m, blk, external_heating=True):

    print(f'\n{"=======> BUILDING BC SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties_feed)
    touch_flow_and_conc(blk.feed)
    blk.product = StateJunction(property_package=m.fs.properties_feed)
    touch_flow_and_conc(blk.product)
    blk.disposal = StateJunction(property_package=m.fs.properties_feed)
    touch_flow_and_conc(blk.disposal)

    blk.recovery_mass = Var(
        initialize=0.5, bounds=(0, 1), units=pyunits.dimensionless, doc="MVC recovery"
    )

    blk.recovery_vol = Var(
        initialize=0.5,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="MVC volumetric recovery",
    )

    # Evaporator
    blk.evaporator = Evaporator(
        property_package_feed=m.fs.properties_feed,
        property_package_vapor=m.fs.properties_vapor,
    )
    # Compressor
    blk.compressor = Compressor(property_package=m.fs.properties_vapor)

    # Condenser
    blk.condenser = Condenser(property_package=m.fs.properties_vapor)

    # Translator SW to Water
    blk.tb_sw_to_water = TranslatorSWtoWater(
        inlet_property_package=m.fs.properties_vapor,
        outlet_property_package=m.fs.properties_feed,
    )

    blk.pump_feed = Pump(property_package=m.fs.properties_feed)
    blk.pump_brine = Pump(property_package=m.fs.properties_feed)
    blk.pump_distillate = Pump(property_package=m.fs.properties_feed)

    blk.separator = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["hx_distillate_cold", "hx_brine_cold"],
        split_basis=SplittingType.totalFlow,
    )

    blk.hx_distillate = HeatExchanger(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.properties_feed, "has_pressure_change": True},
        cold={
            "property_package": m.fs.properties_feed,
            "has_pressure_change": True,
        },
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
    )
    # Set lower bound of approach temperatures
    blk.hx_distillate.delta_temperature_in.setlb(0)
    blk.hx_distillate.delta_temperature_out.setlb(0)
    blk.hx_distillate.area.setlb(5)

    blk.hx_brine = HeatExchanger(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.properties_feed, "has_pressure_change": True},
        cold={
            "property_package": m.fs.properties_feed,
            "has_pressure_change": True,
        },
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
    )
    # Set lower bound of approach temperatures
    blk.hx_brine.delta_temperature_in.setlb(0)
    blk.hx_brine.delta_temperature_out.setlb(0)
    blk.hx_brine.area.setlb(5)

    blk.mixer_feed = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["hx_distillate_cold", "hx_brine_cold"],
    )

    blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump_feed.inlet)
    blk.pump_to_separator = Arc(
        source=blk.pump_feed.outlet, destination=blk.separator.inlet
    )
    blk.sep_dist_cold_to_hx_dist_cold = Arc(
        source=blk.separator.hx_distillate_cold,
        destination=blk.hx_distillate.cold_inlet,
    )
    blk.sep_brine_cold_to_hx_brine_cold = Arc(
        source=blk.separator.hx_brine_cold, destination=blk.hx_brine.cold_inlet
    )
    blk.hx_dist_cold_to_mixer = Arc(
        source=blk.hx_distillate.cold_outlet,
        destination=blk.mixer_feed.hx_distillate_cold,
    )
    blk.hx_brine_cold_to_mixer = Arc(
        source=blk.hx_brine.cold_outlet, destination=blk.mixer_feed.hx_brine_cold
    )
    blk.mixer_feed_to_evaporator = Arc(
        source=blk.mixer_feed.outlet, destination=blk.evaporator.inlet_feed
    )
    blk.evaporator_to_compressor = Arc(
        source=blk.evaporator.outlet_vapor, destination=blk.compressor.inlet
    )
    blk.compressor_to_condenser = Arc(
        source=blk.compressor.outlet, destination=blk.condenser.inlet
    )
    blk.evaporator_to_brine_pump = Arc(
        source=blk.evaporator.outlet_brine, destination=blk.pump_brine.inlet
    )

    blk.brine_pump_to_hx_brine_hot = Arc(
        source=blk.pump_brine.outlet, destination=blk.hx_brine.hot_inlet
    )
    blk.hx_brine_hot_to_disposal = Arc(
        source=blk.hx_brine.hot_outlet, destination=blk.disposal.inlet
    )
    blk.condenser_to_translator = Arc(
        source=blk.condenser.outlet, destination=blk.tb_sw_to_water.inlet
    )
    blk.translated_to_dist_pump = Arc(
        source=blk.tb_sw_to_water.outlet, destination=blk.pump_distillate.inlet
    )

    blk.dist_pump_to_hx_dist_hot = Arc(
        source=blk.pump_distillate.outlet, destination=blk.hx_distillate.hot_inlet
    )
    blk.hx_dist_hot_to_product = Arc(
        source=blk.hx_distillate.hot_outlet, destination=blk.product.inlet
    )

    blk.eq_recovery_mass = Constraint(
        expr=blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
        == blk.recovery_mass
        * (
            blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            + blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
        )
    )

    blk.eq_recovery_vol = Constraint(
        expr=blk.product.properties[0].flow_vol_phase["Liq"]
        == blk.recovery_vol * (blk.feed.properties[0].flow_vol_phase["Liq"])
    )

    blk.eq_separator_split_frac = Constraint(
        expr=blk.separator.split_fraction[0, "hx_distillate_cold"] == blk.recovery_mass
    )

    # blk.hx_area_constr = Constraint(expr=blk.hx_distillate.area >= blk.hx_brine.area)

    if external_heating:
        add_external_heating(blk)

    TransformationFactory("network.expand_arcs").apply_to(blk)

    blk.evaporator.connect_to_condenser(blk.condenser)
    _log.info("MVC flowsheet built")


def set_bc_operating_conditions(
    blk,
    recovery=0.5,
    inlet_brine_temp_guess=50,  # degC
    outlet_brine_temp=70,  # degC
    steam_temp_ub=75,
    compressor_temp_out_lb=65,
    compressor_temp_out_ub=180,
    flow_mass_tds_dist=0.01,
    brine_pump_deltaP=1600,
    dist_pump_deltaP=3700,
    U=3e3,
    **kwargs,
):
    """
    Generic initial point for MVC system.
    """

    blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"].fix(0.1)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(40)
    # Feed pump
    blk.pump_feed.efficiency_pump[0].fix(0.8)
    blk.pump_feed.control_volume.deltaP[0].fix(7e3)

    # Separator
    blk.separator.split_fraction[0, "hx_distillate_cold"] = 0.5

    # Distillate HX
    blk.hx_distillate.overall_heat_transfer_coefficient[0].fix(2e3)
    blk.hx_distillate.area.fix(125)
    blk.hx_distillate.cold.deltaP[0].fix(7e3)
    blk.hx_distillate.hot.deltaP[0].fix(7e3)

    # Brine HX
    blk.hx_brine.overall_heat_transfer_coefficient[0].fix(2e3)
    blk.hx_brine.area.fix(115)
    blk.hx_brine.cold.deltaP[0].fix(7e3)
    blk.hx_brine.hot.deltaP[0].fix(7e3)

    # Evaporator
    blk.evaporator.properties_vapor[0].temperature.setub(steam_temp_ub + 273.15)
    blk.evaporator.inlet_feed.temperature[0] = (
        inlet_brine_temp_guess + 273.15
    )  # provide guess
    blk.evaporator.outlet_brine.temperature[0].fix(outlet_brine_temp + 273.15)
    blk.evaporator.U.fix(U)  # W/K-m^2
    blk.evaporator.area.setub(1e4)  # m^2
    # blk.evaporator.area.set_value(4275)  # m^2

    # Compressor
    blk.compressor.control_volume.properties_out[0].temperature.setlb(
        compressor_temp_out_lb + 273.15
    )
    blk.compressor.control_volume.properties_out[0].temperature.setub(
        compressor_temp_out_ub + 273.15
    )
    blk.compressor.control_volume.properties_out[0].temperature.set_value(
        (compressor_temp_out_lb + compressor_temp_out_ub) / 2 + 273.15
    )
    blk.compressor.pressure_ratio.fix(1.6)
    blk.compressor.efficiency.fix(0.8)

    # Brine pump
    blk.pump_brine.efficiency_pump[0].fix(0.8)
    blk.pump_brine.control_volume.deltaP[0].fix(brine_pump_deltaP)

    # Distillate pump
    blk.pump_distillate.efficiency_pump[0].fix(0.8)
    blk.pump_distillate.control_volume.deltaP[0].fix(dist_pump_deltaP)

    # Fix 0 TDS
    # blk.tb_sw_to_water.properties_out[0].flow_mass_phase_comp["Liq", "TDS"].fix(1e-5)

    blk.tb_sw_to_water.properties_out[0].flow_mass_phase_comp["Liq", "TDS"].fix(
        flow_mass_tds_dist
    )

    print("BC DOF after setting operating conditions: ", degrees_of_freedom(blk))


def set_system_operating_conditions(m):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): m.Qin,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.Cin,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + m.feed_temp,
        },
        hold_state=True,
    )


def set_system_scaling(m):

    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (value(pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s)) * 1000),
        index=("Liq", "H2O"),
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (value(pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s))),
        index=("Liq", "TDS"),
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    set_bc_scaling(m.fs.bc)

    iscale.calculate_scaling_factors(m)


def set_bc_scaling(blk):

    iscale.set_scaling_factor(blk.external_heating, 1e-6)

    # MVC FEED

    # MVC PRODUCT

    # MVC DISPOSAL

    # PUMPS
    iscale.set_scaling_factor(blk.pump_feed.control_volume.work, 1e-3)
    iscale.set_scaling_factor(blk.pump_brine.control_volume.work, 1e-3)
    iscale.set_scaling_factor(blk.pump_distillate.control_volume.work, 1e-3)

    # DISTILLATE HX
    iscale.set_scaling_factor(blk.hx_distillate.hot.heat, 1e-3)
    iscale.set_scaling_factor(blk.hx_distillate.cold.heat, 1e-3)
    iscale.set_scaling_factor(blk.hx_distillate.overall_heat_transfer_coefficient, 1e-3)

    iscale.set_scaling_factor(blk.hx_distillate.area, 1e-1)
    iscale.constraint_scaling_transform(
        blk.hx_distillate.cold_side.pressure_balance[0], 1e-5
    )
    iscale.constraint_scaling_transform(
        blk.hx_distillate.hot_side.pressure_balance[0], 1e-5
    )

    # BRINE HX
    iscale.set_scaling_factor(blk.hx_brine.hot.heat, 1e-3)
    iscale.set_scaling_factor(blk.hx_brine.cold.heat, 1e-3)
    iscale.set_scaling_factor(blk.hx_brine.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(blk.hx_brine.area, 1e-1)
    iscale.constraint_scaling_transform(
        blk.hx_brine.cold_side.pressure_balance[0], 1e-5
    )
    iscale.constraint_scaling_transform(blk.hx_brine.hot_side.pressure_balance[0], 1e-5)

    # EVAPORATOR
    iscale.set_scaling_factor(blk.evaporator.area, 1e-3)
    iscale.set_scaling_factor(blk.evaporator.U, 1e-3)
    iscale.set_scaling_factor(blk.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(blk.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(blk.evaporator.lmtd, 1e-1)

    # COMPRESSOR
    iscale.set_scaling_factor(blk.compressor.control_volume.work, 1e-6)

    # CONDENSER
    iscale.set_scaling_factor(blk.condenser.control_volume.heat, 1e-5)

    iscale.calculate_scaling_factors(blk)


def init_system(m, blk, **kwargs):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_mvc)

    init_bc(blk, **kwargs)

    m.fs.product.initialize()
    propagate_state(m.fs.mvc_to_product)

    m.fs.disposal.initialize()
    propagate_state(m.fs.mvc_to_disposal)


def init_bc(
    # m,
    blk,
    feed_props=None,
    delta_temperature_in=10,
    delta_temperature_out=None,
    solver=None,
    **kwargs,
):
    """
    Initialization routine for generic BC setup.
    To be used with external heating.
    """

    m = blk.model()

    if solver is None:
        solver = get_solver()
    solver.options["halt_on_ampl_error"] = "yes"
    optarg = solver.options

    if feed_props is None:
        feed_props = m.fs.feed.properties[0]

    blk.recovery_mass.fix(0.5)

    blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]

    blk.feed.properties[0].temperature.fix(value(feed_props.temperature))
    blk.feed.properties[0].pressure.fix(value(feed_props.pressure))
    # blk.feed.properties[0].flow_vol_phase["Liq"].fix(value(feed_props.flow_vol_phase["Liq"]))
    # blk.feed.properties[0].conc_mass_phase_comp["Liq", "TDS"].fix(value(feed_props.conc_mass_phase_comp["Liq", "TDS"]))

    solver.solve(blk.feed)

    _log.info(f"{blk.name} feed initialization complete.")

    # Propagate vapor flow rate based on given recovery
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"] = (
        blk.recovery_mass
        * (
            blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            + blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
        )
    )
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Liq", "H2O"] = 0

    # Propagate brine salinity and flow rate
    blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"] = (
        blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
        / (1 - blk.recovery_mass)
    )
    blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "H2O"] = (
        1 - blk.evaporator.properties_brine[0].mass_frac_phase_comp["Liq", "TDS"].value
    )
    blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"] = (
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"] = (
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        - blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
    )

    # Initialize feed pump
    propagate_state(blk.feed_to_pump)
    blk.pump_feed.initialize(
        optarg=optarg,
    )
    _log.info(f"{blk.name} feed pump initialization complete.")

    # Initialize separator
    propagate_state(blk.pump_to_separator)
    # Touch property for initialization
    blk.separator.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.split_fraction[0, "hx_distillate_cold"].fix(blk.recovery_mass.value)
    blk.separator.mixed_state.initialize(
        optarg=optarg,
    )

    # Touch properties for initialization
    blk.separator.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.initialize(
        optarg=optarg,
    )
    blk.separator.split_fraction[0, "hx_distillate_cold"].unfix()
    _log.info(f"{blk.name} separator initialization complete.")

    # Initialize distillate heat exchanger
    propagate_state(blk.sep_dist_cold_to_hx_dist_cold)
    blk.hx_distillate.cold_outlet.temperature[0] = (
        blk.evaporator.inlet_feed.temperature[0].value
    )
    blk.hx_distillate.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[
        0
    ].value
    blk.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value
    )
    blk.hx_distillate.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
    blk.hx_distillate.hot_inlet.temperature[0] = (
        blk.evaporator.outlet_brine.temperature[0].value
    )
    blk.hx_distillate.hot_inlet.pressure[0] = 101325
    blk.hx_distillate.initialize()
    _log.info(f"{blk.name} Distillate HX initialization complete.")

    # Initialize brine heat exchanger
    propagate_state(blk.sep_brine_cold_to_hx_brine_cold)
    blk.hx_brine.cold_outlet.temperature[0] = blk.evaporator.inlet_feed.temperature[
        0
    ].value
    blk.hx_brine.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[0].value
    blk.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    blk.hx_brine.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = (
        blk.evaporator.properties_brine[0].flow_mass_phase_comp["Liq", "TDS"]
    )
    blk.hx_brine.hot_inlet.temperature[0] = blk.evaporator.outlet_brine.temperature[
        0
    ].value
    blk.hx_brine.hot_inlet.pressure[0] = 101325
    blk.hx_brine.initialize()
    _log.info(f"{blk.name} Brine HX initialization complete.")

    # Initialize mixer
    propagate_state(blk.hx_dist_cold_to_mixer)
    propagate_state(blk.hx_brine_cold_to_mixer)
    blk.mixer_feed.initialize()
    blk.mixer_feed.pressure_equality_constraints[0, 2].deactivate()
    _log.info(f"{blk.name} Mixer initialization complete.")

    # Initialize evaporator
    propagate_state(blk.mixer_feed_to_evaporator)
    blk.external_heating.fix()
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
    # fixes and unfixes those values
    blk.evaporator.initialize(
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
    )
    blk.external_heating.unfix()
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    # Initialize compressor
    propagate_state(blk.evaporator_to_compressor)
    blk.compressor.initialize()
    _log.info(f"{blk.name} Compressor initialization complete.")

    # Initialize condenser
    try:
        propagate_state(blk.compressor_to_condenser)
        blk.condenser.initialize(
            heat=-blk.evaporator.heat_transfer.value,
        )
        _log.info(f"{blk.name} Condenser initialization complete.")
    except:
        print_infeasible_constraints(blk.condenser)
        assert False

    # Initialize brine pump
    propagate_state(blk.evaporator_to_brine_pump)
    blk.pump_brine.initialize(
        optarg=optarg,
    )
    _log.info(f"{blk.name} Brine Pump initialization complete.")
    # propagate_state(blk.brine_pump_to_hx_brine_hot)

    # Initialize distillate pump
    propagate_state(blk.condenser_to_translator)  # to translator block
    propagate_state(blk.translated_to_dist_pump)  # from translator block to pump
    blk.pump_distillate.control_volume.properties_in[0].temperature = (
        blk.condenser.control_volume.properties_out[0].temperature.value
    )
    blk.pump_distillate.control_volume.properties_in[0].pressure = (
        blk.condenser.control_volume.properties_out[0].pressure.value
    )
    blk.pump_distillate.initialize(
        optarg=optarg,
    )
    _log.info(f"{blk.name} Distillate Pump initialization complete.")
    # propagate_state(blk.dist_pump_to_hx_dist_hot)

    # Propagate brine state
    propagate_state(blk.hx_brine_hot_to_disposal)
    propagate_state(blk.hx_dist_hot_to_product)

    run_sequential_decomposition(
        blk,
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
        **kwargs,
    )

    blk.product.initialize()
    _log.info(f"{blk.name} Product initialization complete.")
    blk.disposal.initialize()
    _log.info(f"{blk.name} Disposal initialization complete.")

    # m.fs.costing.initialize()
    results = solver.solve(blk)
    _log.info(f"MVC solve termination {results.solver.termination_condition}")
    assert_optimal_termination(results)

    blk.pump_brine.control_volume.deltaP[0].unfix()
    blk.disposal.properties[0].pressure.fix(101325)
    blk.disposal.properties[0].temperature.setub(360)

    print(f"\n~~~~~BC FIRST SOLVE~~~~")
    blk.obj = Objective(expr=blk.external_heating)
    results = solver.solve(blk)
    assert_optimal_termination(results)
    _log.info(f"BC FIRST solve termination {results.solver.termination_condition}")

    blk.external_heating.fix(0)
    blk.del_component(blk.obj)
    blk.evaporator.area.unfix()
    blk.evaporator.outlet_brine.temperature[0].unfix()
    blk.compressor.pressure_ratio.unfix()
    blk.hx_distillate.area.unfix()
    blk.hx_brine.area.unfix()

    print(f"BC dof = {degrees_of_freedom(blk)}")
    print(f"model dof = {degrees_of_freedom(m)}")
    results = solver.solve(blk)
    assert_optimal_termination(results)
    print(f"BC SECOND solve termination {results.solver.termination_condition}")
    _log.info(f"BC SECOND solve termination {results.solver.termination_condition}")

    mvc_feed_state_vars = blk.feed.properties[0].define_port_members()
    feed_state_vars = feed_props.define_port_members()
    blk.feed.properties[0].mass_frac_phase_comp.unfix()
    for k, v in mvc_feed_state_vars.items():
        if v.is_indexed():
            for i, vv in v.items():
                vv.fix(value(feed_state_vars[k][i]))
        else:
            v.fix(value(feed_state_vars[k]))
    try:
        print(f"\n~~~~~BC THIRD SOLVE~~~~")
        print(f"BC dof = {degrees_of_freedom(blk)}")
        print(f"model dof = {degrees_of_freedom(m)}")
        blk.feed.initialize()
        # m.fs.costing.initialize()
        results = solver.solve(blk)
        assert_optimal_termination(results)
        _log.info(f"BC THIRD solve termination {results.solver.termination_condition}")
    except:
        print(f"BC THIRD SOLVE FAILED!!\n")
        pass

    mvc_feed_state_vars = blk.feed.properties[0].define_port_members()
    feed_state_vars = feed_props.define_port_members()
    for k, v in mvc_feed_state_vars.items():
        if v.is_indexed():
            for i, vv in v.items():
                vv.unfix()
        else:
            v.unfix()

    blk.feed.initialize()
    blk.feed.properties[0].temperature.unfix()
    blk.feed.properties[0].pressure.unfix()

    print(f"Initialization done, BC dof = {degrees_of_freedom(blk)}")
    print(f"Initialization done, model dof = {degrees_of_freedom(m)}")


def run_sequential_decomposition(
    blk,
    delta_temperature_in=60,
    delta_temperature_out=None,
    tear_solver="cbc",
    iterlim=5,
):

    def func_initialize(unit):
        if unit.local_name in ["feed", "product", "disposal"]:
            pass
        elif unit.local_name == "condenser":
            unit.initialize(
                heat=-unit.flowsheet().evaporator.heat_transfer.value,
                optarg=solver.options,
            )
        elif unit.local_name == "evaporator":
            unit.flowsheet().external_heating.fix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
            unit.initialize(
                delta_temperature_in=delta_temperature_in,
                delta_temperature_out=delta_temperature_out,
            )
            unit.flowsheet().external_heating.unfix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
        elif unit.local_name == "separator":
            unit.split_fraction[0, "hx_distillate_cold"].fix(
                unit.flowsheet().recovery_mass.value
            )
            unit.initialize()
            unit.split_fraction[0, "hx_distillate_cold"].unfix()
        elif unit.local_name == "mixer_feed":
            unit.initialize()
            unit.pressure_equality_constraints[0, 2].deactivate()
        else:
            unit.initialize()

    seq = SequentialDecomposition(tear_solver=tear_solver)
    seq.options.log_info = True
    seq.options.iterLim = iterlim
    seq.run(blk, func_initialize)


def add_external_heating(blk):
    """
    Add additional heat source so infeasible designs can be used
    as initial guess prior to optimization
    """

    blk.external_heating = Var(
        initialize=0,
        units=pyunits.watt,
        bounds=(0, None),
        doc="External heating for evaporator",
    )

    blk.evaporator.eq_energy_balance.deactivate()
    blk.evaporator.eq_energy_balance_with_external_heat = Constraint(
        expr=blk.evaporator.heat_transfer
        + blk.external_heating
        + blk.evaporator.properties_feed[0].enth_flow
        == blk.evaporator.properties_brine[0].enth_flow
        + blk.evaporator.properties_vapor[0].enth_flow_phase["Vap"]
    )
    iscale.set_scaling_factor(blk.external_heating, 1e-6)


def add_bc_costing(blk, flowsheet_costing_block=None):

    if flowsheet_costing_block is None:
        m = blk.model()
        flowsheet_costing_block = m.fs.costing

    blk.evaporator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.compressor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )

    blk.pump_feed.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.pump_distillate.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.pump_brine.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.hx_distillate.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.hx_brine.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )
    blk.mixer_feed.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )

    blk.costing = Block()
    blk.costing.capital_cost = Expression(
        expr=sum(
            pyunits.convert(cap, to_units=flowsheet_costing_block.base_currency)
            for cap in [
                blk.evaporator.costing.capital_cost
                + blk.compressor.costing.capital_cost
                + blk.pump_feed.costing.capital_cost
                + blk.pump_distillate.costing.capital_cost
                + blk.pump_brine.costing.capital_cost
                + blk.hx_brine.costing.capital_cost
                + blk.hx_distillate.costing.capital_cost
                + blk.mixer_feed.costing.capital_cost
            ]
        )
    )

    # m.fs.costing.TIC.fix(2)
    # m.fs.costing.electricity_cost = 0.1  # 0.15
    flowsheet_costing_block.heat_exchanger.material_factor_cost.fix(5)
    flowsheet_costing_block.evaporator.material_factor_cost.fix(5)
    flowsheet_costing_block.compressor.unit_cost.fix(1 * 7364)


def set_up_optimization(m, blk):
    """
    Add objective and unfix design variables
    """
    if hasattr(blk, "LCOW_obj"):
        blk.del_component(blk.LCOW_obj)
    blk.LCOW_obj = Objective(expr=m.fs.costing.LCOW)
    blk.external_heating.fix(0)
    blk.evaporator.area.unfix()
    blk.evaporator.outlet_brine.temperature[0].unfix()
    blk.compressor.pressure_ratio.unfix()
    blk.hx_distillate.area.unfix()
    blk.hx_brine.area.unfix()


def display_bc_flow_table(blk, w=25):
    title = "BC System Flow Table"
    side = int(((5 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(
        f'{"Unit":<{w}s}{"Mass Flow Water (kg/s)":<{w}s}{"Pressure (bar)":<{w}s}{"Mass Flow NaCl (kg/s)":<{w}s}{"Conc. (g/L)":<{w}s}'
    )
    print(f"{'-' * (5 * w)}")
    print(
        f'{"Feed":<{w}s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<{w}.2f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )
    print(
        f'{"Product":<{w}s}{blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(blk.product.properties[0].pressure, to_units=pyunits.bar)():<{w}.2f}{blk.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )
    print(
        f'{"Disposal":<{w}s}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(blk.disposal.properties[0].pressure, to_units=pyunits.bar)():<{w}.2f}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )

    print("\n\n")


def calculate_cost_sf(cost):
    print(cost.name, cost.value)
    if cost.value in [0, None]:
        sf = 1e-2
    else:
        sf = 10 ** -(math.log10(abs(cost.value)))
    iscale.set_scaling_factor(cost, sf)


def scale_bc_costs(m, blk):
    calculate_cost_sf(blk.hx_distillate.costing.capital_cost)
    calculate_cost_sf(blk.hx_brine.costing.capital_cost)
    calculate_cost_sf(blk.mixer_feed.costing.capital_cost)
    calculate_cost_sf(blk.evaporator.costing.capital_cost)
    calculate_cost_sf(blk.compressor.costing.capital_cost)
    # calculate_cost_sf(m.fs.costing.aggregate_capital_cost)
    # calculate_cost_sf(m.fs.costing.aggregate_flow_costs["electricity"])
    # calculate_cost_sf(m.fs.costing.total_capital_cost)
    # calculate_cost_sf(m.fs.costing.total_operating_cost)

    iscale.calculate_scaling_factors(m)

    print("Scaled costs")


def print_bc_stream_flows(blk, w=30):

    flow_in = pyunits.convert(
        blk.feed.properties[0.0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.min,
    )
    flow_out = pyunits.convert(
        blk.product.properties[0.0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.min,
    )
    flow_brine = pyunits.convert(
        blk.disposal.properties[0.0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.min,
    )
    conc_in = pyunits.convert(
        blk.feed.properties[0.0].conc_mass_phase_comp["Liq", "TDS"],
        to_units=pyunits.g / pyunits.L,
    )
    conc_out = pyunits.convert(
        blk.product.properties[0.0].conc_mass_phase_comp["Liq", "TDS"],
        to_units=pyunits.g / pyunits.L,
    )
    conc_brine = pyunits.convert(
        blk.disposal.properties[0.0].conc_mass_phase_comp["Liq", "TDS"],
        to_units=pyunits.g / pyunits.L,
    )
    print(f"{'Flow In':<{w}s}{value(flow_in):<{w}.3f}{'gpm':<{w}s}")
    print(f"{'Flow Distillate':<{w}s}{value(flow_out):<{w}.3f}{'gpm':<{w}s}")
    print(f"{'Flow Brine':<{w}s}{value(flow_brine):<{w}.3f}{'gpm':<{w}s}")
    print(f"{'Conc TDS In':<{w}s}{value(conc_in):<{w}.3f}{'g/L':<{w}s}")
    print(f"{'Conc TDS Distillate':<{w}s}{value(conc_out):<{w}.3f}{'g/L':<{w}s}")
    print(f"{'Conc TDS Brine':<{w}s}{value(conc_brine):<{w}.3f}{'g/L':<{w}s}")


def report_pump(blk, w=35):

    pump_power_watt = value(pyunits.convert(blk.work_mechanical[0], to_units=pyunits.W))
    # pump_flow_in = pyunits.convert(
    #     blk.control_volume.properties_in[0].flow_vol_phase["Liq"],
    #     to_units=pyunits.gallon / pyunits.min,
    # )
    pump_flow_out = pyunits.convert(
        blk.control_volume.properties_out[0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.min,
    )
    pressure_in = pyunits.convert(
        blk.control_volume.properties_in[0].pressure, to_units=pyunits.psi
    )
    pressure_out = pyunits.convert(
        blk.control_volume.properties_out[0].pressure, to_units=pyunits.psi
    )
    # print(f'{"Pump Inlet Flow":<{w}s}{value(pump_flow_in):<{w}.2f}{"gpm":<{w}s}')
    print(f'{"Pump Flow":<{w}s}{value(pump_flow_out):<{w}.2f}{"gpm":<{w}s}')

    print(f'{"Pump Power":<{w}s}{value(pump_power_watt):<{w}.2f}{"W":<{w}s}')
    print(
        f'{"DeltaP":<{w}s}{value(blk.deltaP[0]):<{w}.2f}{f"{pyunits.get_units(blk.deltaP[0])}":<{w}s}'
    )
    print(
        f'{"Pressure In":<{w}s}{value(pressure_in):<{w}.2f}{f"{pyunits.get_units(pressure_in)}":<{w}s}'
    )
    print(
        f'{"Pressure Out":<{w}s}{value(pressure_out):<{w}.2f}{f"{pyunits.get_units(pressure_out)}":<{w}s}'
    )
    print(
        f'{"Temp. In":<{w}s}{value(blk.control_volume.properties_out[0].temperature):<{w}.2f}{f"{pyunits.get_units(blk.control_volume.properties_out[0].temperature)}":<{w}s}'
    )
    print(
        f'{"Temp. Out":<{w}s}{value(blk.control_volume.properties_out[0].temperature):<{w}.2f}{f"{pyunits.get_units(blk.control_volume.properties_out[0].temperature)}":<{w}s}'
    )
    temp_out_F = (
        blk.control_volume.properties_out[0].temperature.value - 273.15
    ) * 9 / 5 + 32
    print(f'{"Temp. Out (F)":<{w}s}{temp_out_F:<{w}.2f}{"F":<{w}s}')
    print(
        f'{"Pressure Ratio":<{w}s}{value(blk.ratioP[0]):<{w}.2f}{f"{pyunits.get_units(blk.ratioP[0])}":<{w}s}'
    )


def report_bc(blk, w=35):
    title = "MVC System Flow Table"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print(f'\n{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print_bc_stream_flows(blk, w=w)

    dist_temp_F = (blk.product.properties[0].temperature.value - 273.15) * 9 / 5 + 32
    print(
        f"{'Distillate temperature':<{w}s}{blk.product.properties[0].temperature.value:<{w}.3f}{f'{pyunits.get_units(blk.product.properties[0].temperature)}':<{w}s}"
    )
    print(f"{'Distillate temperature (F)':<{w}s}{dist_temp_F:<{w}.3f}{'F':<{w}s}")
    print(f'{"Recovery (vol)":<{w}s}{blk.recovery_vol.value:<{w}.3f}{"gal/gal":<{w}s}')
    print(f'{"Recovery (mass)":<{w}s}{blk.recovery_mass.value:<{w}.3f}{"kg/kg":<{w}s}')

    if hasattr(blk, "pump_feed"):
        print(f"\nFeed Pump {'.' * w}")
        report_pump(blk.pump_feed, w=w)
    if hasattr(blk, "compressor"):
        print(f"\nCompressor {'.' * w}")
        temp_out_F = (
            blk.compressor.control_volume.properties_out[0].temperature.value - 273.15
        ) * 9 / 5 + 32
        pressure_in_psi = value(
            pyunits.convert(
                blk.compressor.control_volume.properties_in[0].pressure,
                to_units=pyunits.psi,
            )
        )
        pressure_out_psi = value(
            pyunits.convert(
                blk.compressor.control_volume.properties_out[0].pressure,
                to_units=pyunits.psi,
            )
        )
        print(
            f'{"Inlet Pressure":<{w}s}{blk.compressor.control_volume.properties_in[0].pressure.value:<{w}.3f}{f"{pyunits.get_units(blk.compressor.control_volume.properties_in[0].pressure)}":<{w}s}'
        )
        print(f'{"Inlet Pressure (psi)":<{w}s}{pressure_in_psi:<{w}.3f}{"psi":<{w}s}')
        print(
            f'{"Outlet Pressure":<{w}s}{blk.compressor.control_volume.properties_out[0].pressure.value:<{w}.3f}{f"{pyunits.get_units(blk.compressor.control_volume.properties_out[0].pressure)}":<{w}s}'
        )
        print(f'{"Outlet Pressure (psi)":<{w}s}{pressure_out_psi:<{w}.3f}{"psi":<{w}s}')

        print(
            f'{"Vapor Temp.":<{w}s}{blk.compressor.control_volume.properties_out[0].temperature.value:<{w}.3f}{f"{pyunits.get_units(blk.compressor.control_volume.properties_out[0].temperature)}":<{w}s}'
        )
        print(f'{"Vapor Temp. (F)":<{w}s}{temp_out_F:<{w}.3f}{"F":<{w}s}')
        print(
            f'{"Vapor Pressure":<{w}s}{blk.compressor.control_volume.properties_out[0].pressure.value:<{w}.3f}{f"{pyunits.get_units(blk.compressor.control_volume.properties_out[0].pressure)}":<{w}s}'
        )
        print(
            f'{"Pressure Ratio":<{w}s}{blk.compressor.pressure_ratio.value:<{w}.3f}{f"{pyunits.get_units(blk.compressor.pressure_ratio)}":<{w}s}'
        )

    if hasattr(blk, "condenser"):
        cond_pressure_psi = value(
            pyunits.convert(
                blk.condenser.control_volume.properties_out[0].pressure,
                to_units=pyunits.psi,
            )
        )
        print(f"\nCondenser {'.' * w}")
        print(
            f'{"Vapor Temp.":<{w}s}{blk.condenser.control_volume.properties_out[0].temperature.value:<{w}.3f}{f"{pyunits.get_units(blk.condenser.control_volume.properties_out[0].temperature)}":<{w}s}'
        )
        print(
            f'{"Vapor Pressure":<{w}s}{blk.condenser.control_volume.properties_out[0].pressure.value:<{w}.3f}{f"{pyunits.get_units(blk.condenser.control_volume.properties_out[0].pressure)}":<{w}s}'
        )
        print(f'{"Vapor Pressure (psi)":<{w}s}{cond_pressure_psi:<{w}.3f}{"psi":<{w}s}')
        condensed_vap_temp_F = (
            blk.condenser.control_volume.properties_out[0].temperature.value - 273.15
        ) * 9 / 5 + 32
        print(
            f'{"Condensed Vapor Temp. (F)":<{w}s}{condensed_vap_temp_F:<{w}.3f}{"F":<{w}s}'
        )
    if hasattr(blk, "evaporator"):
        print(f"\nEvaporator {'.' * w}")
        print(
            f'{"Preheated Feed Temp.":<{w}s}{blk.evaporator.properties_feed[0].temperature.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.properties_feed[0].temperature)}":<{w}s}'
        )
        print(
            f'{"Brine Pressure":<{w}s}{blk.evaporator.properties_brine[0].pressure.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.properties_brine[0].pressure)}":<{w}s}'
        )
        print(
            f'{"Brine/Vapor Temp.":<{w}s}{blk.evaporator.properties_brine[0].temperature.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.properties_brine[0].temperature)}":<{w}s}'
        )
        # print(
        #     f'{"Vapor Temp.":<{w}s}{blk.evaporator.properties_vapor[0].temperature.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.properties_vapor[0].temperature)}":<{w}s}'
        # )
        print(
            f'{"Evaporator Area":<{w}s}{blk.evaporator.area.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.area)}":<{w}s}'
        )
        print(
            f'{"LMTD":<{w}s}{blk.evaporator.lmtd.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.lmtd)}":<{w}s}'
        )
        print(
            f'{"Delta T_in":<{w}s}{blk.evaporator.delta_temperature_in.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.delta_temperature_in)}":<{w}s}'
        )
        print(
            f'{"Delta T_out":<{w}s}{blk.evaporator.delta_temperature_out.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.delta_temperature_out)}":<{w}s}'
        )
        U2 = pyunits.convert(
            blk.evaporator.U,
            to_units=pyunits.BTU / (pyunits.hr * pyunits.ft**2 * pyunits.degR),
        )
        print(
            f'{"U":<{w}s}{blk.evaporator.U.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.U)}":<{w}s}'
        )
        print(f'{"U (BTU/hr-ft2-R)":<{w}s}{value(U2):<{w}.2f}{"BTU/hr-ft2-R":<{w}s}')
        # print(
        #     f'{"Duty":<{w}s}{blk.evaporator.heat_transfer.value:<{w}.2f}{f"{pyunits.get_units(blk.evaporator.heat_transfer)}":<{w}s}'
        # )
    if hasattr(blk, "hx_brine"):
        print(f"\nBrine Heat Exchanger {'.' * w}")
        print(
            f'{"Overall HX Coeff":<{w}s}{blk.hx_brine.overall_heat_transfer_coefficient[0].value:<{w}.2f}{f"{pyunits.get_units(blk.hx_brine.overall_heat_transfer_coefficient[0])}":<{w}s}'
        )
        print(
            f'{"HX Brine Area":<{w}s}{blk.hx_brine.area.value:<{w}.2f}{f"{pyunits.get_units(blk.hx_brine.area)}":<{w}s}'
        )
        print(
            f'{"Delta T_in":<{w}s}{blk.hx_brine.delta_temperature_in[0].value:<{w}.2f}{f"{pyunits.get_units(blk.hx_brine.delta_temperature_in[0])}":<{w}s}'
        )
        print(
            f'{"Delta T_out":<{w}s}{blk.hx_brine.delta_temperature_out[0].value:<{w}.2f}{f"{pyunits.get_units(blk.hx_brine.delta_temperature_out[0])}":<{w}s}'
        )
    if hasattr(blk, "hx_distillate"):
        print(f"\nDistillate Heat Exchanger {'.' * w}")
        print(
            f'{"Overall HX Coeff":<{w}s}{blk.hx_distillate.overall_heat_transfer_coefficient[0].value:<{w}.2f}{f"{pyunits.get_units(blk.hx_distillate.overall_heat_transfer_coefficient[0])}":<{w}s}'
        )
        print(
            f'{"HX Distillate Area":<{w}s}{blk.hx_distillate.area.value:<{w}.2f}{f"{pyunits.get_units(blk.hx_distillate.area)}":<{w}s}'
        )
        print(
            f'{"Delta T_in":<{w}s}{blk.hx_distillate.delta_temperature_in[0].value:<{w}.2f}{f"{pyunits.get_units(blk.hx_distillate.delta_temperature_in[0])}":<{w}s}'
        )
        print(
            f'{"Delta T_out":<{w}s}{blk.hx_distillate.delta_temperature_out[0].value:<{w}.2f}{f"{pyunits.get_units(blk.hx_distillate.delta_temperature_out[0])}":<{w}s}'
        )
    if hasattr(blk, "pump_brine"):
        print(f"\nBrine Pump {'.' * w}")
        report_pump(blk.pump_brine, w=w)
    if hasattr(blk, "pump_distillate"):
        print(f"\nDistillate Pump {'.' * w}")
        report_pump(blk.pump_distillate, w=w)

    print("\n\n")


def main(recovery_vol=0.95):

    m = build_system()
    set_system_scaling(m)
    set_system_operating_conditions(m)
    set_bc_operating_conditions(m.fs.bc)
    init_system(m, m.fs.bc)
    m.fs.bc.recovery_vol.fix(recovery_vol)
    m.fs.bc.recovery_mass.unfix()
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)
    print(f"dof = {degrees_of_freedom(m)}")
    _ = solve_bc(m)
    propagate_state(m.fs.bc.dist_pump_to_hx_dist_hot)
    propagate_state(m.fs.bc.brine_pump_to_hx_brine_hot)
    m.fs.bc.mixer_feed.pressure_equality_constraints.activate()
    print(f"dof = {degrees_of_freedom(m)}")
    _ = solve_bc(m)
    print(f"dof = {degrees_of_freedom(m)}")
    report_bc(m.fs.bc)
    m.fs.costing.LCOW.display()
    m.fs.costing.SEC.display()

    return m


if __name__ == "__main__":
    m = main()
    # m = main(recovery_vol=0.5)
