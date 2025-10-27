import pathlib
import math

from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Param,
    Var,
    Constraint,
    Expression,
    Objective,
    check_optimal_termination,
    assert_optimal_termination,
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
    Heater,
)
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.pressure_changer import Pump
from watertap.unit_models.mvc.components import Evaporator
from watertap.unit_models.mvc.components import Compressor
from watertap.unit_models.mvc.components import Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap.costing import WaterTAPCosting

from srp.components.translator_sw_to_water import Translator_SW_to_Water

__author__ = "Kurban Sitterley"

__all__ = [
    "build_mvc",
    "build_and_run_mvc",
    "set_mvc_operating_conditions",
    "set_mvc_scaling",
    "init_mvc",
    "add_mvc_costing",
    "scale_mvc_costs",
    "run_sequential_decomposition",
    "solve_mvc",
    "display_MVC_flow_table",
    "scale_mvc_costs",
    "report_pump",
    "report_MVC",
]

solver = get_solver()

_log = idaeslog.getLogger("SRP")


def build_and_run_mvc(
    recovery=0.5,
    Qin=728 / 30,  # gpm
    Cin=11.408,
    **kwargs,
):
    Qin = Qin * pyunits.gallons / pyunits.minute
    # Cin_wells = 11408 * pyunits.mg / pyunits.L
    Cin = Cin * pyunits.g / pyunits.liter
    feed_temp = 27  # degrees C
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties_feed)
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)() * 1000),
        index=("Liq", "H2O"),
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)()),
        index=("Liq", "TDS"),
    )

    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase
    iscale.calculate_scaling_factors(m)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin,
            ("conc_mass_phase_comp", ("Liq", "TDS")): Cin_wells,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + feed_temp,
        },
        hold_state=True,
    )
    m.fs.BC = FlowsheetBlock()
    build_mvc(m, m.fs.BC)
    add_mvc_costing(m, m.fs.BC)
    m.fs.costing.cost_process()
    set_mvc_operating_conditions(m, m.fs.BC)
    set_mvc_scaling(m, m.fs.BC)
    init_mvc(m, m.fs.BC)

    m.fs.feed_to_BC = Arc(source=m.fs.feed.outlet, destination=m.fs.BC.feed.inlet)
    propagate_state(m.fs.feed_to_BC)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    print(f"dof = {degrees_of_freedom(m)}")
    m.fs.feed.initialize()
    m.fs.BC.recovery_vol.fix(0.9)
    m.fs.BC.recovery_mass.unfix()
    m.fs.BC.feed.properties[0].flow_vol_phase
    m.fs.BC.feed.properties[0].conc_mass_phase_comp
    m.fs.BC.feed.initialize()
    m.fs.BC.product.properties[0].flow_vol_phase
    m.fs.BC.product.properties[0].conc_mass_phase_comp
    m.fs.BC.product.initialize()
    m.fs.BC.disposal.properties[0].flow_vol_phase
    m.fs.BC.disposal.properties[0].conc_mass_phase_comp
    m.fs.BC.disposal.initialize()
    results = solve_mvc(m)
    print_infeasible_constraints(m)

    report_MVC(m.fs.BC)
    print(f"dof = {degrees_of_freedom(m)}")

    return m


def solve_mvc(blk):
    solver = get_solver()
    try:
        results = solver.solve(blk, tee=False)
        print(f"termination MVC {results.solver.termination_condition}")
        assert_optimal_termination(results)
    except:
        results = solver.solve(blk, tee=False)
        print(f"termination MVC {results.solver.termination_condition}")
        assert_optimal_termination(results)

    return results


def build_mvc_system(recovery=0.5, **kwargs):

    m = ConcreteModel()
    m.recovery_mass = recovery
    m.recovery_vol = recovery

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties_feed)
    m.fs.product = Product(property_package=m.fs.properties_feed)
    m.fs.disposal = Product(property_package=m.fs.properties_feed)

    m.fs.MVC = mvc = FlowsheetBlock(dynamic=False)

    build_mvc(m, mvc, **kwargs)

    m.fs.feed_to_mvc = Arc(source=m.fs.feed.outlet, destination=mvc.feed.inlet)

    m.fs.mvc_to_product = Arc(source=mvc.product.outlet, destination=m.fs.product.inlet)

    m.fs.mvc_to_disposal = Arc(
        source=mvc.disposal.outlet, destination=m.fs.disposal.inlet
    )

    add_mvc_costing(m, mvc)
    # print(pyunits.get_units(mvc.evaporator.costing.capital_cost))
    # assert False
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(mvc.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(mvc.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        mvc.product.properties[0].flow_vol, name="SEC"
    )

    m.fs.costing.heat_exchanger.material_factor_cost.fix(5)
    m.fs.costing.evaporator.material_factor_cost.fix(5)
    m.fs.costing.compressor.unit_cost.fix(1 * 7364)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_mvc(m, blk, external_heating=True):

    print(f'\n{"=======> BUILDING MVC SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties_feed)
    blk.product = StateJunction(property_package=m.fs.properties_feed)
    blk.disposal = StateJunction(property_package=m.fs.properties_feed)

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
    blk.tb_sw_to_water = Translator_SW_to_Water(
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

    if external_heating:
        add_external_heating(m, blk)

    TransformationFactory("network.expand_arcs").apply_to(m)

    blk.evaporator.connect_to_condenser(blk.condenser)
    _log.info("MVC flowsheet built")


def set_mvc_operating_conditions(
    m,
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
    print("DOF after setting operating conditions: ", degrees_of_freedom(blk))


def set_system_operating_conditions(m, Qin=1, tds=130, feed_temp=25):

    global flow_in

    Qin = Qin * pyunits.Mgallons / pyunits.day
    tds = tds * pyunits.g / pyunits.liter
    flow_in = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_in,
            ("conc_mass_phase_comp", ("Liq", "TDS")): tds,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + feed_temp,
        },
        hold_state=True,
    )


def set_mvc_scaling(
    m,
    blk,
    properties_feed=None,
    properties_vapor=None,
    calc_blk_scaling_factors=True,
):

    if properties_feed is None:
        properties_feed = m.fs.properties_feed

    if properties_vapor is None:
        properties_vapor = m.fs.properties_vapor

    properties_feed.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )

    iscale.set_scaling_factor(blk.external_heating, 1e-6)

    # MVC FEED
    # set_scaling_factor(
    #     blk.feed.properties[0.0].conc_mass_phase_comp["Liq", "TDS"], 1e-2
    # )
    # set_scaling_factor(blk.hx_distillate.hot_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"], 1e2)

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
    # set_scaling_factor(blk.evaporator.properties_brine[0].pressure, 1e-4)
    iscale.set_scaling_factor(blk.evaporator.area, 1e-3)
    iscale.set_scaling_factor(blk.evaporator.U, 1e-3)
    iscale.set_scaling_factor(blk.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(blk.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(blk.evaporator.lmtd, 1e-1)
    # set_scaling_factor(blk.evaporator.heat_transfer, 1e-7)

    # COMPRESSOR
    iscale.set_scaling_factor(blk.compressor.control_volume.work, 1e-6)
    # set_scaling_factor(
    #     blk.compressor.control_volume.properties_in[0].enth_flow_phase["Liq"], 1
    # )

    # CONDENSER
    iscale.set_scaling_factor(blk.condenser.control_volume.heat, 1e-5)

    # if hasattr(blk.evaporator, "costing"):
    #     scale_mvc_costs(m, blk)

    if calc_blk_scaling_factors:
        iscale.calculate_scaling_factors(blk)

    else:
        iscale.calculate_scaling_factors(m)


def init_system(m, blk, **kwargs):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_mvc)

    init_mvc(m, blk, **kwargs)

    m.fs.product.initialize()
    propagate_state(m.fs.mvc_to_product)

    m.fs.disposal.initialize()
    propagate_state(m.fs.mvc_to_disposal)


def init_mvc(
    m,
    blk,
    feed_props=None,
    delta_temperature_in=10,
    delta_temperature_out=None,
    solver=None,
    **kwargs,
):
    """
    Initialization routine for generic MVC setup.
    To be used with external heating.
    """

    if solver is None:
        solver = get_solver()
    solver.options["halt_on_ampl_error"] = "yes"
    # solver.options["tee"] = True
    optarg = solver.options

    if feed_props is None:
        feed_props = m.fs.feed.properties[0]

    blk.recovery_mass.fix(0.5)

    blk.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]

    blk.feed.properties[0].temperature.fix(value(feed_props.temperature))
    blk.feed.properties[0].pressure.fix(value(feed_props.pressure))

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
    blk.pump_feed.initialize(optarg=optarg, solver="ipopt-watertap")
    _log.info(f"{blk.name} feed pump initialization complete.")

    # Initialize separator
    propagate_state(blk.pump_to_separator)
    # Touch property for initialization
    blk.separator.mixed_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.split_fraction[0, "hx_distillate_cold"].fix(blk.recovery_mass.value)
    blk.separator.mixed_state.initialize(optarg=optarg, solver="ipopt-watertap")

    # Touch properties for initialization
    blk.separator.hx_brine_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.hx_distillate_cold_state[0].mass_frac_phase_comp["Liq", "TDS"]
    blk.separator.initialize(optarg=optarg, solver="ipopt-watertap")
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
    blk.hx_distillate.initialize(solver="ipopt-watertap")
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
    blk.hx_brine.initialize(solver="ipopt-watertap")
    _log.info(f"{blk.name} Brine HX initialization complete.")

    # Initialize mixer
    propagate_state(blk.hx_dist_cold_to_mixer)
    propagate_state(blk.hx_brine_cold_to_mixer)
    blk.mixer_feed.initialize(solver="ipopt-watertap")
    blk.mixer_feed.pressure_equality_constraints[0, 2].deactivate()
    _log.info(f"{blk.name} Mixer initialization complete.")
    # print(f"\ndof @ 1 = {degrees_of_freedom(blk)}\n")

    # Initialize evaporator
    propagate_state(blk.mixer_feed_to_evaporator)
    blk.external_heating.fix()
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
    # fixes and unfixes those values

    blk.evaporator.initialize(
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
        solver="ipopt-watertap",
    )
    blk.external_heating.unfix()
    blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

    # Initialize compressor
    propagate_state(blk.evaporator_to_compressor)
    blk.compressor.initialize(solver="ipopt-watertap")
    _log.info(f"{blk.name} Compressor initialization complete.")

    # Initialize condenser
    try:
        propagate_state(blk.compressor_to_condenser)
        blk.condenser.initialize(
            heat=-blk.evaporator.heat_transfer.value, solver="ipopt-watertap"
        )
        _log.info(f"{blk.name} Condenser initialization complete.")
    except:
        print_infeasible_constraints(blk.condenser)
        assert False

    # Initialize brine pump
    propagate_state(blk.evaporator_to_brine_pump)
    blk.pump_brine.initialize(optarg=optarg, solver="ipopt-watertap")
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
    blk.pump_distillate.initialize(optarg=optarg, solver="ipopt-watertap")
    _log.info(f"{blk.name} Distillate Pump initialization complete.")
    # propagate_state(blk.dist_pump_to_hx_dist_hot)

    # Propagate brine state
    propagate_state(blk.hx_brine_hot_to_disposal)
    propagate_state(blk.hx_dist_hot_to_product)

    # print(f"\ndof @ 2 = {degrees_of_freedom(blk)}\n")

    run_sequential_decomposition(
        m,
        blk,
        delta_temperature_in=delta_temperature_in,
        delta_temperature_out=delta_temperature_out,
        **kwargs,
    )

    blk.product.initialize()
    _log.info(f"{blk.name} Product initialization complete.")
    blk.disposal.initialize()
    _log.info(f"{blk.name} Disposal initialization complete.")

    m.fs.costing.initialize()
    results = solver.solve(blk, tee=False)
    print(f"MVC solve termination {results.solver.termination_condition}")
    _log.info(f"MVC solve termination {results.solver.termination_condition}")
    assert_optimal_termination(results)

    print(f"blk dof at end of init = {degrees_of_freedom(blk)}")

    blk.pump_brine.control_volume.deltaP[0].unfix()
    blk.disposal.properties[0].pressure.fix(101325)
    blk.disposal.properties[0].temperature.setub(360)

    print(f"\n~~~~~FIRST SOLVE~~~~")
    blk.obj = Objective(expr=blk.external_heating)
    print(f"blk dof = {degrees_of_freedom(blk)}")
    print(f"model dof = {degrees_of_freedom(m)}")
    results = solver.solve(blk, tee=False)
    assert_optimal_termination(results)
    print(f"MVC FIRST solve termination {results.solver.termination_condition}")
    _log.info(f"MVC FIRST solve termination {results.solver.termination_condition}")

    blk.external_heating.fix(0)
    del blk.obj
    blk.external_heating.fix(0)
    blk.evaporator.area.unfix()
    blk.evaporator.outlet_brine.temperature[0].unfix()
    blk.compressor.pressure_ratio.unfix()
    blk.hx_distillate.area.unfix()
    blk.hx_brine.area.unfix()

    print(f"\n~~~~~SECOND SOLVE~~~~")

    print(f"blk dof = {degrees_of_freedom(blk)}")
    print(f"model dof = {degrees_of_freedom(m)}")
    results = solver.solve(blk, tee=False)
    assert_optimal_termination(results)
    print(f"MVC SECOND solve termination {results.solver.termination_condition}")
    _log.info(f"MVC SECOND solve termination {results.solver.termination_condition}")

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
        print(f"\n~~~~~THIRD SOLVE~~~~")
        print(f"blk dof = {degrees_of_freedom(blk)}")
        print(f"model dof = {degrees_of_freedom(m)}")
        blk.feed.initialize()
        m.fs.costing.initialize()
        results = solver.solve(blk, tee=False)
        # results = solve_mvc(blk)
        assert_optimal_termination(results)
        print(f"MVC THIRD solve termination {results.solver.termination_condition}")
        _log.info(f"MVC THIRD solve termination {results.solver.termination_condition}")
    except:
        print(f"MVC THIRD SOLVE FAILED!!\n")
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

    print(f"Initialization done, blk dof = {degrees_of_freedom(blk)}")
    print(f"Initialization done, model dof = {degrees_of_freedom(m)}")


def run_sequential_decomposition(
    m,
    blk,
    delta_temperature_in=60,
    delta_temperature_out=None,
    tear_solver="cbc",
    iterlim=5,
):

    def func_initialize(unit):
        # print(unit.local_name)
        # print(f"dof = {degrees_of_freedom(unit)}\n")
        # if unit.local_name == "feed":
        if unit.local_name in ["feed", "product", "disposal"]:
            pass
        elif unit.local_name == "condenser":
            unit.initialize(
                heat=-unit.flowsheet().evaporator.heat_transfer.value,
                optarg=solver.options,
                solver="ipopt-watertap",
            )
        elif unit.local_name == "evaporator":
            unit.flowsheet().external_heating.fix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
            unit.initialize(
                delta_temperature_in=delta_temperature_in,
                delta_temperature_out=delta_temperature_out,
                solver="ipopt-watertap",
            )
            unit.flowsheet().external_heating.unfix()
            unit.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()
        elif unit.local_name == "separator":
            unit.split_fraction[0, "hx_distillate_cold"].fix(
                unit.flowsheet().recovery_mass.value
            )
            unit.initialize(solver="ipopt-watertap")
            unit.split_fraction[0, "hx_distillate_cold"].unfix()
        elif unit.local_name == "mixer_feed":
            unit.initialize(solver="ipopt-watertap")
            unit.pressure_equality_constraints[0, 2].deactivate()
        # elif unit.local_name == "hx_distillate":
        #     unit.cold_outlet.temperature[0] = blk.evaporator.inlet_feed.temperature[
        #         0
        #     ].value
        #     unit.cold_outlet.pressure[0] = blk.evaporator.inlet_feed.pressure[0].value
        #     unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] = (
        #         blk.evaporator.properties_vapor[0]
        #         .flow_mass_phase_comp["Vap", "H2O"]
        #         .value
        #     )
        #     unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"] = 1e-4
        #     unit.hot_inlet.temperature[0] = blk.evaporator.outlet_brine.temperature[
        #         0
        #     ].value
        #     unit.hot_inlet.pressure[0] = 101325
        #     unit.initialize(solver="ipopt-watertap")
        else:
            unit.initialize(solver="ipopt-watertap")

    seq = SequentialDecomposition(tear_solver=tear_solver)
    seq.options.log_info = True
    seq.options.iterLim = iterlim
    seq.run(blk, func_initialize)


def add_external_heating(m, blk):
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


def add_mvc_costing(m, blk, flowsheet_costing_block=None):

    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing

    flowsheet_costing_block.base_currency = pyunits.USD_2023

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
            pyunits.convert(cap, to_units=pyunits.USD_2023)
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


def display_MVC_flow_table(blk, w=25):
    title = "MVC System Flow Table"
    side = int(((5 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(
        f'{"Unit":<{w}s}{"Mass Flow Water (kg/s)":<{w}s}{"Pressure (bar)":<{w}s}{"Mass Flow NaCl (kg/s)":<{w}s}{"Conc. (g/L)":<{w}s}'
    )
    print(f"{'-' * (5 * w)}")
    print(
        f'{"Feed":<{w}s}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{value(pyunits.convert(blk.feed.properties[0.0].pressure, to_units=pyunits.bar)):<{w}.1f}{blk.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )
    print(
        f'{"Product":<{w}s}{blk.product.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(blk.product.properties[0].pressure, to_units=pyunits.bar)():<{w}.1f}{blk.product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )
    print(
        f'{"Disposal":<{w}s}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value:<{w}.3f}{pyunits.convert(blk.disposal.properties[0].pressure, to_units=pyunits.bar)():<{w}.1f}{blk.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value:<{w}.3e}{blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<{w}.3f}'
    )

    print("\n\n")


def calculate_cost_sf(cost):
    print(cost.name, cost.value)
    if cost.value in [0, None]:
        sf = 1e-2
    else:
        sf = 10 ** -(math.log10(abs(cost.value)))
    iscale.set_scaling_factor(cost, sf)


def scale_mvc_costs(m, blk):
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


def print_MVC_stream_flows(blk, w=30):
    # sb_in = model.fs.BCs.find_component(f"to_{blk.name}_state")
    # conc_in = sb_in[0].conc_mass_phase_comp["Liq", "TDS"]
    # print(f"{'Flow in ':<{w}s}{value(pyunits.convert(flow_in, to_units=pyunits.gallon / pyunits.min)):<{w}.3f}{'gpm':<{w}s}")
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
    # print(f'{"Pump Inlet Flow":<{w}s}{value(pump_flow_in):<{w}.1f}{"gpm":<{w}s}')
    print(f'{"Pump Flow":<{w}s}{value(pump_flow_out):<{w}.1f}{"gpm":<{w}s}')

    print(f'{"Pump Power":<{w}s}{value(pump_power_watt):<{w}.1f}{"W":<{w}s}')
    print(
        f'{"DeltaP":<{w}s}{value(blk.deltaP[0]):<{w}.1f}{f"{pyunits.get_units(blk.deltaP[0])}":<{w}s}'
    )
    print(
        f'{"Pressure In":<{w}s}{value(pressure_in):<{w}.1f}{f"{pyunits.get_units(pressure_in)}":<{w}s}'
    )
    print(
        f'{"Pressure Out":<{w}s}{value(pressure_out):<{w}.1f}{f"{pyunits.get_units(pressure_out)}":<{w}s}'
    )
    print(
        f'{"Temp. In":<{w}s}{value(blk.control_volume.properties_out[0].temperature):<{w}.1f}{f"{pyunits.get_units(blk.control_volume.properties_out[0].temperature)}":<{w}s}'
    )
    print(
        f'{"Temp. Out":<{w}s}{value(blk.control_volume.properties_out[0].temperature):<{w}.1f}{f"{pyunits.get_units(blk.control_volume.properties_out[0].temperature)}":<{w}s}'
    )
    temp_out_F = (
        blk.control_volume.properties_out[0].temperature.value - 273.15
    ) * 9 / 5 + 32
    print(f'{"Temp. Out (F)":<{w}s}{temp_out_F:<{w}.1f}{"F":<{w}s}')
    print(
        f'{"Pressure Ratio":<{w}s}{value(blk.ratioP[0]):<{w}.1f}{f"{pyunits.get_units(blk.ratioP[0])}":<{w}s}'
    )


def report_MVC(blk, w=35):
    title = "MVC System Flow Table"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print(f'\n{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print_MVC_stream_flows(blk, w=w)

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
            f'{"Preheated Feed Temp.":<{w}s}{blk.evaporator.properties_feed[0].temperature.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.properties_feed[0].temperature)}":<{w}s}'
        )
        print(
            f'{"Brine Pressure":<{w}s}{blk.evaporator.properties_brine[0].pressure.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.properties_brine[0].pressure)}":<{w}s}'
        )
        print(
            f'{"Brine/Vapor Temp.":<{w}s}{blk.evaporator.properties_brine[0].temperature.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.properties_brine[0].temperature)}":<{w}s}'
        )
        # print(
        #     f'{"Vapor Temp.":<{w}s}{blk.evaporator.properties_vapor[0].temperature.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.properties_vapor[0].temperature)}":<{w}s}'
        # )
        print(
            f'{"Evaporator Area":<{w}s}{blk.evaporator.area.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.area)}":<{w}s}'
        )
        print(
            f'{"LMTD":<{w}s}{blk.evaporator.lmtd.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.lmtd)}":<{w}s}'
        )
        print(
            f'{"Delta T_in":<{w}s}{blk.evaporator.delta_temperature_in.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.delta_temperature_in)}":<{w}s}'
        )
        print(
            f'{"Delta T_out":<{w}s}{blk.evaporator.delta_temperature_out.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.delta_temperature_out)}":<{w}s}'
        )
        U2 = pyunits.convert(
            blk.evaporator.U,
            to_units=pyunits.BTU / (pyunits.hr * pyunits.ft**2 * pyunits.degR),
        )
        print(
            f'{"U":<{w}s}{blk.evaporator.U.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.U)}":<{w}s}'
        )
        print(f'{"U (BTU/hr-ft2-R)":<{w}s}{value(U2):<{w}.1f}{"BTU/hr-ft2-R":<{w}s}')
        # print(
        #     f'{"Duty":<{w}s}{blk.evaporator.heat_transfer.value:<{w}.1f}{f"{pyunits.get_units(blk.evaporator.heat_transfer)}":<{w}s}'
        # )
    if hasattr(blk, "hx_brine"):
        print(f"\nBrine Heat Exchanger {'.' * w}")
        print(
            f'{"Overall HX Coeff":<{w}s}{blk.hx_brine.overall_heat_transfer_coefficient[0].value:<{w}.1f}{f"{pyunits.get_units(blk.hx_brine.overall_heat_transfer_coefficient[0])}":<{w}s}'
        )
        print(
            f'{"HX Brine Area":<{w}s}{blk.hx_brine.area.value:<{w}.1f}{f"{pyunits.get_units(blk.hx_brine.area)}":<{w}s}'
        )
        print(
            f'{"Delta T_in":<{w}s}{blk.hx_brine.delta_temperature_in[0].value:<{w}.1f}{f"{pyunits.get_units(blk.hx_brine.delta_temperature_in[0])}":<{w}s}'
        )
        print(
            f'{"Delta T_out":<{w}s}{blk.hx_brine.delta_temperature_out[0].value:<{w}.1f}{f"{pyunits.get_units(blk.hx_brine.delta_temperature_out[0])}":<{w}s}'
        )
    if hasattr(blk, "hx_distillate"):
        print(f"\nDistillate Heat Exchanger {'.' * w}")
        print(
            f'{"Overall HX Coeff":<{w}s}{blk.hx_distillate.overall_heat_transfer_coefficient[0].value:<{w}.1f}{f"{pyunits.get_units(blk.hx_distillate.overall_heat_transfer_coefficient[0])}":<{w}s}'
        )
        print(
            f'{"HX Distillate Area":<{w}s}{blk.hx_distillate.area.value:<{w}.1f}{f"{pyunits.get_units(blk.hx_distillate.area)}":<{w}s}'
        )
        print(
            f'{"Delta T_in":<{w}s}{blk.hx_distillate.delta_temperature_in[0].value:<{w}.1f}{f"{pyunits.get_units(blk.hx_distillate.delta_temperature_in[0])}":<{w}s}'
        )
        print(
            f'{"Delta T_out":<{w}s}{blk.hx_distillate.delta_temperature_out[0].value:<{w}.1f}{f"{pyunits.get_units(blk.hx_distillate.delta_temperature_out[0])}":<{w}s}'
        )
    if hasattr(blk, "pump_brine"):
        print(f"\nBrine Pump {'.' * w}")
        report_pump(blk.pump_brine, w=w)
    if hasattr(blk, "pump_distillate"):
        print(f"\nDistillate Pump {'.' * w}")
        report_pump(blk.pump_distillate, w=w)

    print("\n\n")


if __name__ == "__main__":
    Qin = 350 * pyunits.gallons / pyunits.minute
    Cin_wells = 11408 * pyunits.mg / pyunits.L
    feed_temp = 27  # degrees C
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties_feed)
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)() * 1000),
        index=("Liq", "H2O"),
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)()),
        index=("Liq", "TDS"),
    )

    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase
    iscale.calculate_scaling_factors(m)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin,
            ("conc_mass_phase_comp", ("Liq", "TDS")): Cin_wells,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + feed_temp,
        },
        hold_state=True,
    )
    m.fs.BC = FlowsheetBlock()
    build_mvc(m, m.fs.BC)
    add_mvc_costing(m, m.fs.BC)
    m.fs.costing.cost_process()
    set_mvc_operating_conditions(m, m.fs.BC)
    set_mvc_scaling(m, m.fs.BC)
    init_mvc(m, m.fs.BC)

    m.fs.feed_to_BC = Arc(source=m.fs.feed.outlet, destination=m.fs.BC.feed.inlet)
    propagate_state(m.fs.feed_to_BC)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    print(f"dof = {degrees_of_freedom(m)}")
    m.fs.feed.initialize()
    m.fs.BC.recovery_vol.fix(0.95)
    m.fs.BC.recovery_mass.unfix()
    m.fs.BC.feed.properties[0].flow_vol_phase
    m.fs.BC.feed.properties[0].conc_mass_phase_comp
    m.fs.BC.feed.initialize()
    m.fs.BC.product.properties[0].flow_vol_phase
    m.fs.BC.product.properties[0].conc_mass_phase_comp
    m.fs.BC.product.initialize()
    m.fs.BC.disposal.properties[0].flow_vol_phase
    m.fs.BC.disposal.properties[0].conc_mass_phase_comp
    m.fs.BC.disposal.initialize()
    m.fs.BC.hx_area_constr = Constraint(
        expr=m.fs.BC.hx_distillate.area <= m.fs.BC.hx_brine.area
    )
    results = solve_mvc(m)
    print_infeasible_constraints(m)

    report_MVC(m.fs.BC)
    print(f"dof = {degrees_of_freedom(m)}")

    # m.fs.BC.compressor.pressure_ratio.fix(1.6)
    # results = solve_mvc(m)
    # print_infeasible_constraints(m)

    # report_MVC(m.fs.BC)
    # print(f'dof = {degrees_of_freedom(m)}')

    # m.fs.BC.hx_distillate.area.fix()
    # results = solve_mvc(m)
    # print_infeasible_constraints(m)

    # report_MVC(m.fs.BC)
    # print(f'dof = {degrees_of_freedom(m)}')

    # m.fs.BC.hx_brine.area.fix()
    # results = solve_mvc(m)
    # print_infeasible_constraints(m)

    # report_MVC(m.fs.BC)
    # print(f'dof = {degrees_of_freedom(m)}')

    # m.fs.BC.evaporator.area.fix()
    # results = solve_mvc(m)
    # print_infeasible_constraints(m)

    # report_MVC(m.fs.BC)
    # print(f'dof = {degrees_of_freedom(m)}')

    m.fs.costing.add_LCOW(m.fs.BC.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.BC.product.properties[0].flow_vol_phase["Liq"], name="SEC"
    )
    m.fs.obj = Objective(expr=m.fs.costing.LCOW)
    results = solve_mvc(m)
    print_infeasible_constraints(m)

    report_MVC(m.fs.BC)
    print(f"dof = {degrees_of_freedom(m)}")
    m.fs.costing.SEC.display()
    m.fs.costing.LCOW.display()

    # m.fs.BC.evaporator.outlet_brine.temperature[0].unfix()
    # m.fs.BC.compressor.pressure_ratio.unfix()
    # print(f'dof = {degrees_of_freedom(m)}')
    # for x in [0.91, 0.92, 0.93, 0.94, 0.95,]:
    #     m.fs.BC.recovery_vol.fix(x)
    #     run_sequential_decomposition(m, m.fs.BC, iterlim=25)
    #     try:
    #         results = solve_mvc(m)
    #     except:
    #         print_infeasible_constraints(m)

    #     # report_MVC(m.fs.BC)
    #     print(f'dof = {degrees_of_freedom(m)}')
    # m.fs.BC.recovery_vol.fix(0.91)
    # try:
    #     results = solve_mvc(m)
    # except:
    #     print_infeasible_constraints(m)

    # report_MVC(m.fs.BC)

    # print(f"LCOW = {value(m.fs.costing.LCOW):.3f} $/m3")

    # m.fs.BC.pump_brine.display()
    # # m = build_and_run_mvc(recovery=0.48, Qin=4.9, tds=117.416)
    # U2 = pyunits.convert(m.fs.BC.evaporator.U, to_units=pyunits.BTU / (pyunits.hr * pyunits.ft**2 * pyunits.degR))
    # hx2 = pyunits.convert(m.fs.BC.evaporator.heat_transfer, to_units=pyunits.BTU / pyunits.hr)
    # print(f"U = {U2():.1f} BTU/(hr*ft^2*degR)")
    # print(f"HX = {hx2():.1f} BTU/hr")
    # m.fs.BC.compressor.control_volume.properties_out[0].display()

    # m = build_and_run_mvc(recovery=0.9, Qin=0.424920699462, tds=11.408)
    # blk = m.fs.MVC
    # recov = value(blk.product.properties[0].flow_vol_phase["Liq"])
    # m.fs.feed.properties[0].flow_vol_phase.display()
    # m.fs.feed.properties[0].conc_mass_phase_comp.display()
    # m.fs.product.properties[0].conc_mass_phase_comp.display()
    # m.fs.disposal.properties[0].conc_mass_phase_comp.display()
    # blk.recovery_vol.display()
    # blk.evaporator.costing.display()
    # blk.compressor.costing.display()

    # blk.pump_feed.costing.display()
    # blk.pump_distillate.costing.display()
    # blk.pump_brine.costing.display()
    # blk.hx_distillate.costing.display()
    # blk.hx_brine.costing.display()
    # blk.mixer_feed.costing.display()
    # blk.costing.capital_cost.display()
    # m.fs.MVC.recovery_vol.display()
    # m.fs.costing.total_capital_cost.display()
    # m.fs.costing.total_operating_cost.display()
    # m.fs.product.properties[0].flow_vol_phase.display()
    # m.fs.MVC.product.properties[0].flow_vol_phase.display()
    # m.fs.costing.LCOW.display()


## USED TO RE-INIT MVC @ THIRD SOLVE
# blk.condenser.initialize(heat=-blk.evaporator.heat_transfer.value)

# blk.external_heating.fix()
# blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].fix()
# blk.evaporator.initialize(
#     delta_temperature_in=delta_temperature_in,
#     delta_temperature_out=delta_temperature_out,
#     solver="ipopt-watertap",
# )
# blk.external_heating.unfix()
# blk.evaporator.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].unfix()

# blk.separator.split_fraction[0, "hx_distillate_cold"].fix(
#     blk.recovery_mass.value
# )
# blk.separator.initialize(solver="ipopt-watertap")
# blk.separator.split_fraction[0, "hx_distillate_cold"].unfix()

# blk.mixer_feed.initialize(solver="ipopt-watertap")
# blk.mixer_feed.pressure_equality_constraints[0, 2].deactivate()


# m = ConcreteModel()
# m.fs = FlowsheetBlock()

# m.fs.properties_feed = SeawaterParameterBlock()
# m.fs.properties_vapor = SteamParameterBlock()

# m.fs.wells = Feed(property_package=m.fs.properties_feed)

# m.fs.raw_water_tank = Separator(
#     property_package=m.fs.properties_feed,
#     outlet_list=["to_cooling_tower", "to_service_and_fire", "to_uf"],
#     split_basis=SplittingType.componentFlow,
# )
# m.fs.raw_water_tank.split_fraction[0, "to_uf", "H2O"].fix(114.2 / 11343)
# m.fs.raw_water_tank.split_fraction[0, "to_uf", "TDS"].fix(114.2 / 11343)
# m.fs.raw_water_tank.split_fraction[0, "to_service_and_fire", "H2O"].fix(17 / 11343)
# m.fs.raw_water_tank.split_fraction[0, "to_service_and_fire", "TDS"].fix(17 / 11343)
# # m.fs.raw_water_tank.to_service_and_fire_state[0].conc_mass_phase_comp
# # m.fs.raw_water_tank.to_uf_state[0].conc_mass_phase_comp

# m.fs.service_and_fire = Product(property_package=m.fs.properties_feed)

# m.fs.uf = Separator(
#     property_package=m.fs.properties_feed,
#     outlet_list=["to_ro", "to_ro_containment"],
#     split_basis=SplittingType.componentFlow,
# )
# m.fs.uf.split_fraction[0, "to_ro", "H2O"].fix(91.3 / 114.2)
# m.fs.uf.split_fraction[0, "to_ro", "TDS"].fix(91.3 / 114.2)

# m.fs.ro = Separator(
#     property_package=m.fs.properties_feed,
#     outlet_list=["to_ro_containment", "to_ro_permeate"],
#     split_basis=SplittingType.componentFlow,
# )
# m.fs.ro.split_fraction[0, "to_ro_permeate", "H2O"].fix(56.2 / 91.3)
# m.fs.ro.split_fraction[0, "to_ro_permeate", "TDS"].fix(0.01)
# m.fs.ro.to_ro_permeate_state[0].conc_mass_phase_comp
# m.fs.ro.to_ro_containment_state[0].conc_mass_phase_comp

# m.fs.permeate = Product(property_package=m.fs.properties_feed)

# m.fs.cooling_tower = Separator(
#     property_package=m.fs.properties_feed,
#     outlet_list=["to_evaporation", "to_ww_surge_tank"],
#     split_basis=SplittingType.componentFlow,
# )
# m.fs.cooling_tower.split_fraction[0, "to_evaporation", "H2O"].fix(9178.9 / 10522.4)
# m.fs.cooling_tower.split_fraction[0, "to_evaporation", "TDS"].fix(0)
# m.fs.cooling_tower.to_ww_surge_tank_state[0].conc_mass_phase_comp

# m.fs.evaporation = Product(property_package=m.fs.properties_feed)

# m.fs.ww_surge_tank = Separator(
#     property_package=m.fs.properties_feed,
#     outlet_list=["to_ro_reject_tank", "to_bc_feed_tank"],
#     split_basis=SplittingType.componentFlow,
# )
# m.fs.ww_surge_tank.split_fraction[0, "to_bc_feed_tank", "H2O"].fix(544.5 / 1353.1)
# m.fs.ww_surge_tank.split_fraction[0, "to_bc_feed_tank", "TDS"].fix(544.5 / 1353.1)

# m.fs.bc_feed_tank = Mixer(
#     property_package=m.fs.properties_feed,
#     momentum_mixing_type=MomentumMixingType.none,
#     inlet_list=["from_ro_reject_tank", "from_ww_surge_tank"],
# )
# m.fs.bc_feed_tank.outlet.pressure[0].fix(101325)

# m.fs.ro_reject_tank = Separator(
#     property_package=m.fs.properties_feed,
#     outlet_list=["to_bc_feed_tank", "to_ro_containment", "to_conc_waste"],
#     split_basis=SplittingType.componentFlow,
# )
# m.fs.ro_reject_tank.split_fraction[0, "to_ro_containment", "H2O"].fix(115.2 / 305.1)
# m.fs.ro_reject_tank.split_fraction[0, "to_ro_containment", "TDS"].fix(115.2 / 305.1)
# m.fs.ro_reject_tank.split_fraction[0, "to_conc_waste", "H2O"].fix(6.4 / 305.1)
# m.fs.ro_reject_tank.split_fraction[0, "to_conc_waste", "TDS"].fix(6.4 / 305.1)
# m.fs.ro_reject_tank.to_ro_containment_state[0].conc_mass_phase_comp
# m.fs.ro_reject_tank.to_conc_waste_state[0].conc_mass_phase_comp
# m.fs.ro_reject_tank.to_bc_feed_tank_state[0].conc_mass_phase_comp
# m.fs.ro_reject_tank.mixed_state[0].conc_mass_phase_comp

# m.fs.ro_containment = Mixer(
#     property_package=m.fs.properties_feed,
#     momentum_mixing_type=MomentumMixingType.none,
#     inlet_list=["from_uf", "from_ro_reject_tank", "from_ro"],
# )
# m.fs.ro_containment.outlet.pressure[0].fix(101325)

# # m.fs.ro_reject_tank.to_conc_waste_state[0].conc_mass_phase_comp
# # m.fs.ro_reject_tank.to_ro_containment_state[0].conc_mass_phase_comp

# # m.fs.conc_waste = Mixer(
# #     property_package=m.fs.properties_feed,
# #     momentum_mixing_type=MomentumMixingType.equality,
# #     # inlet_list=["from_bc", "from_ro_reject_tank"],
# #     inlet_list=["from_ro_reject_tank"],
# # )
# m.fs.conc_waste = StateJunction(property_package=m.fs.properties_feed)

# m.fs.bc = Product(property_package=m.fs.properties_feed)

# m.fs.evaporation_ponds = Mixer(
#     property_package=m.fs.properties_feed,
#     momentum_mixing_type=MomentumMixingType.none,
#     inlet_list=["from_conc_waste", "from_ro_containment"],
# )
# m.fs.evaporation_ponds.outlet.pressure[0].fix(101325)


# m.fs.evaporation_2 = Product(property_package=m.fs.properties_feed)
# # =======================================================================
# # Connections
# # =======================================================================
# m.fs.wells_to_raw_water_tank = Arc(
#     source=m.fs.wells.outlet, destination=m.fs.raw_water_tank.inlet
# )

# m.fs.raw_water_tank_to_cooling_tower = Arc(
#     source=m.fs.raw_water_tank.to_cooling_tower, destination=m.fs.cooling_tower.inlet
# )
# m.fs.raw_water_tank_to_service_and_fire = Arc(
#     source=m.fs.raw_water_tank.to_service_and_fire,
#     destination=m.fs.service_and_fire.inlet,
# )
# m.fs.raw_water_tank_to_uf = Arc(
#     source=m.fs.raw_water_tank.to_uf, destination=m.fs.uf.inlet
# )

# m.fs.cooling_tower_to_evaporation = Arc(
#     source=m.fs.cooling_tower.to_evaporation, destination=m.fs.evaporation.inlet
# )
# m.fs.cooling_tower_to_ww_surge_tank = Arc(
#     source=m.fs.cooling_tower.to_ww_surge_tank, destination=m.fs.ww_surge_tank.inlet
# )

# m.fs.ww_surge_tank_to_ro_reject_tank = Arc(
#     source=m.fs.ww_surge_tank.to_ro_reject_tank, destination=m.fs.ro_reject_tank.inlet
# )
# m.fs.ww_surge_tank_to_bc_feed_tank = Arc(
#     source=m.fs.ww_surge_tank.to_bc_feed_tank,
#     destination=m.fs.bc_feed_tank.from_ww_surge_tank,
# )

# m.fs.ro_reject_tank_to_bc_feed_tank = Arc(
#     source=m.fs.ro_reject_tank.to_bc_feed_tank,
#     destination=m.fs.bc_feed_tank.from_ro_reject_tank,
# )
# m.fs.ro_reject_tank_to_ro_containment = Arc(
#     source=m.fs.ro_reject_tank.to_ro_containment,
#     destination=m.fs.ro_containment.from_ro_reject_tank,
# )
# m.fs.ro_reject_tank_to_conc_waste = Arc(
#     source=m.fs.ro_reject_tank.to_conc_waste,
#     # destination=m.fs.conc_waste.from_ro_reject_tank,
#     destination=m.fs.conc_waste.inlet,
# )

# m.fs.uf_to_ro_containment = Arc(
#     source=m.fs.uf.to_ro_containment, destination=m.fs.ro_containment.from_uf
# )
# m.fs.uf_to_ro = Arc(source=m.fs.uf.to_ro, destination=m.fs.ro.inlet)

# m.fs.ro_to_ro_containment = Arc(
#     source=m.fs.ro.to_ro_containment, destination=m.fs.ro_containment.from_ro
# )
# m.fs.ro_to_ro_permeate = Arc(
#     source=m.fs.ro.to_ro_permeate, destination=m.fs.permeate.inlet
# )

# m.fs.conc_waste_to_evaporation_ponds = Arc(
#     source=m.fs.conc_waste.outlet, destination=m.fs.evaporation_ponds.from_conc_waste
# )

# m.fs.ro_containment_to_evaporation_ponds = Arc(
#     source=m.fs.ro_containment.outlet, destination=m.fs.evaporation_ponds.from_ro_containment
# )

# m.fs.evaporation_ponds_to_evaporation_2 = Arc(
#     source=m.fs.evaporation_ponds.outlet, destination=m.fs.evaporation_2.inlet
# )

# m.fs.bc_feed_tank_to_bc = Arc(
#     source=m.fs.bc_feed_tank.outlet, destination=m.fs.bc.inlet
# )
# TransformationFactory("network.expand_arcs").apply_to(m)
# # =======================================================================
# # Initialization
# # =======================================================================

# Qin = 11343 * pyunits.gallons / pyunits.minute
# Cin_wells = 1467 * pyunits.mg / pyunits.L
# feed_temp = 27  # degrees C

# # =======================================================================
# # Scaling
# # =======================================================================

# m.fs.properties_feed.set_default_scaling(
#     "flow_mass_phase_comp",
#     1 / (pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)() * 1000),
#     index=("Liq", "H2O"),
# )
# m.fs.properties_feed.set_default_scaling(
#     "flow_mass_phase_comp",
#     1 / (pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)() * 1000),
#     index=("Liq", "TDS"),
# )
# m.fs.properties_feed.set_default_scaling(
#     "flow_mass_phase_comp",
#     1 / (pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)() * 1000),
#     index=("Liq", "H2O"),
# )
# m.fs.properties_feed.set_default_scaling(
#     "flow_vol_phase",
#     1 / pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)(),
#     index=("Liq",),
# )
# m.fs.wells.properties[0].conc_mass_phase_comp
# m.fs.wells.properties[0].flow_vol_phase
# iscale.calculate_scaling_factors(m)

# m.fs.wells.properties.calculate_state(
#     var_args={
#         ("flow_vol_phase", ("Liq")): Qin,
#         ("conc_mass_phase_comp", ("Liq", "TDS")): Cin_wells,
#         ("pressure", None): 101325,
#         ("temperature", None): 273.15 + feed_temp,
#     },
#     hold_state=True,
# )

# print(f"dof = {degrees_of_freedom(m)}")
# m.fs.wells.initialize()
# propagate_state(m.fs.wells_to_raw_water_tank)

# m.fs.raw_water_tank.initialize()
# propagate_state(m.fs.raw_water_tank_to_cooling_tower)
# propagate_state(m.fs.raw_water_tank_to_service_and_fire)
# propagate_state(m.fs.raw_water_tank_to_uf)

# m.fs.cooling_tower.initialize()
# propagate_state(m.fs.cooling_tower_to_evaporation)
# propagate_state(m.fs.cooling_tower_to_ww_surge_tank)
# m.fs.evaporation.initialize()

# m.fs.service_and_fire.initialize()

# m.fs.uf.initialize()
# propagate_state(m.fs.uf_to_ro_containment)
# propagate_state(m.fs.uf_to_ro)

# m.fs.ro.initialize()
# propagate_state(m.fs.ro_to_ro_containment)
# propagate_state(m.fs.ro_to_ro_permeate)

# m.fs.permeate.initialize()

# m.fs.ww_surge_tank.initialize()
# propagate_state(m.fs.ww_surge_tank_to_ro_reject_tank)
# propagate_state(m.fs.ww_surge_tank_to_bc_feed_tank)

# m.fs.ro_reject_tank.initialize()
# propagate_state(m.fs.ro_reject_tank_to_ro_containment)
# propagate_state(m.fs.ro_reject_tank_to_bc_feed_tank)
# propagate_state(m.fs.ro_reject_tank_to_conc_waste)

# m.fs.bc_feed_tank.mixed_state[0].flow_vol_phase
# m.fs.bc_feed_tank.mixed_state[0].conc_mass_phase_comp
# m.fs.bc_feed_tank.initialize()
# propagate_state(m.fs.bc_feed_tank_to_bc)

# m.fs.bc.initialize()

# m.fs.conc_waste.initialize()
# propagate_state(m.fs.conc_waste_to_evaporation_ponds)
# m.fs.ro_containment.initialize()
# propagate_state(m.fs.conc_waste_to_evaporation_ponds)
# propagate_state(m.fs.ro_containment_to_evaporation_ponds)

# m.fs.evaporation_ponds.initialize()
# propagate_state(m.fs.evaporation_ponds_to_evaporation_2)

# m.fs.evaporation_2.initialize()

# print(f"dof = {degrees_of_freedom(m)}")

# # results = solver.solve(m, tee=True)
