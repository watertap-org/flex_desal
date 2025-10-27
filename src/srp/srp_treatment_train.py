import pandas as pd
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Param,
    Objective,
    Block,
    Var,
    TransformationFactory,
    units as pyunits,
    check_optimal_termination,
    assert_optimal_termination,
)
from pyomo.network import Arc, SequentialDecomposition

import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Separator, Mixer, Product, StateJunction
from idaes.models.unit_models.translator import Translator
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale

from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.unit_models.pressure_changer import Pump
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
from watertap.costing import WaterTAPCosting
import math
import matplotlib.pyplot as plt
from watertap.core.util.model_diagnostics.infeasible import (
    print_infeasible_constraints,
    print_infeasible_bounds,
    print_variables_close_to_bounds,
    print_constraints_close_to_bounds,
    print_close_to_bounds,
)

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)

import srp.components.brine_concentrator as bc

solver = get_solver()


def build_primary_fs(
    RO_pump_pressure=350,
    BCs=["BC_C", "BC_A", "BC_B"],
    Qin=11343,
    Cin=1467,
    feed_temp=27,
):

    Qin = Qin * pyunits.gallons / pyunits.minute
    Cin = Cin * pyunits.mg / pyunits.L
    # feed_temp = 27  # degrees C

    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.costing = WaterTAPCosting()

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    m.BCs = BCs
    # ------------------------ WELLS ------------------------
    m.fs.wells = Feed(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.wells)

    m.fs.raw_water_tank_mix = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["from_permeate", "from_wells"],
    )
    m.fs.raw_water_tank_mix.outlet.pressure[0].fix(101325)
    touch_flow_and_conc(m.fs.raw_water_tank_mix)

    # ------------------------ RAW WATER TANK ------------------------
    m.fs.raw_water_tank = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["to_cooling_tower", "to_service_and_fire"],
        split_basis=SplittingType.componentFlow,
    )
    m.fs.raw_water_tank.split_fraction[0, "to_cooling_tower", "H2O"].fix(
        10522.4 / 11343
    )
    m.fs.raw_water_tank.split_fraction[0, "to_cooling_tower", "TDS"].fix(
        10522.4 / 11343
    )
    touch_flow_and_conc(m.fs.raw_water_tank)

    # ------------------------ SERVICE AND FIRE ------------------------
    m.fs.service_and_fire = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.service_and_fire)

    # ------------------------ UF SYSTEM ------------------------
    m.fs.uf = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["to_ro_pump", "to_ro_containment"],
        split_basis=SplittingType.componentFlow,
    )
    m.fs.uf.split_fraction[0, "to_ro_pump", "H2O"].fix(91.3 / 114.2)
    m.fs.uf.split_fraction[0, "to_ro_pump", "TDS"].fix(91.3 / 114.2)
    touch_flow_and_conc(m.fs.uf)

    # ------------------------ RO PUMP ------------------------
    m.fs.ro_pump = Pump(
        property_package=m.fs.properties_feed,
    )
    m.fs.ro_pump.control_volume.properties_out[0].pressure.fix(
        RO_pump_pressure * pyunits.psi
    )
    # m.fs.ro_pump.costing = UnitModelCostingBlock(
    #     flowsheet_costing_block=m.fs.costing,
    # )
    m.fs.ro_pump.efficiency_pump.fix(0.8)

    m.fs.ro_pump.op_pressure_slope = Param(
        initialize=15.7,
        mutable=True,
        units=pyunits.psi / (pyunits.kg / pyunits.m**3),
        doc="Slope of pump outlet pressure vs salinity",
    )
    m.fs.ro_pump.op_pressure_intercept = Param(
        initialize=96.47,
        mutable=True,
        units=pyunits.psi,
        doc="Intercept of pump outlet pressure vs salinity",
    )
    m.fs.ro_pump.pressure_constr = Constraint(
        expr=m.fs.ro_pump.control_volume.properties_out[0].pressure
        == pyunits.convert(
            m.fs.ro_pump.op_pressure_slope
            * m.fs.ro_pump.control_volume.properties_in[0].conc_mass_phase_comp[
                "Liq", "TDS"
            ]
            + m.fs.ro_pump.op_pressure_intercept,
            to_units=pyunits.Pa,
        )
    )
    m.fs.ro_pump.pressure_constr.deactivate()

    touch_flow_and_conc(m.fs.ro_pump)

    # ------------------------ RO ------------------------
    m.fs.feed_osm_pressure = Param(
        initialize=9.76,
        units=pyunits.bar,
        mutable=True,
        doc="Feed osmotic pressure",  # estimated from vant Hoff equation for NaCl @ 25C
    )
    m.fs.ro = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["to_ro_containment", "to_ro_permeate"],
        split_basis=SplittingType.componentFlow,
    )
    # m.fs.ro.split_fraction[0, "to_ro_permeate", "H2O"].fix(56.2 / 91.3)
    m.fs.ro.split_fraction[0, "to_ro_permeate", "TDS"].fix(0.025)
    touch_flow_and_conc(m.fs.ro)

    m.fs.ro.recov_base = Param(initialize=0.1370946, mutable=True)
    m.fs.ro.recov_exp = Param(initialize=0.58863718, mutable=True)
    m.fs.ro_recovery_constr = Constraint(
        expr=m.fs.ro.split_fraction[0, "to_ro_permeate", "H2O"]
        == m.fs.ro.recov_base
        * (
            pyunits.convert(
                m.fs.ro_pump.control_volume.properties_out[0].pressure,
                to_units=pyunits.bar,
            )
            - m.fs.feed_osm_pressure
        )
        ** m.fs.ro.recov_exp
    )
    m.fs.ro.split_fraction[0, "to_ro_permeate", "H2O"].setlb(0.5)
    m.fs.ro.split_fraction[0, "to_ro_permeate", "H2O"].setub(0.8)

    # ------------------------ PERMEATE ------------------------
    m.fs.permeate = StateJunction(property_package=m.fs.properties_feed)
    m.fs.permeate = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["to_demin", "to_raw_water_tank_mix"],
        split_basis=SplittingType.componentFlow,
    )
    # m.fs.permeate.split_fraction[0, "to_demin", "H2O"].fix(1.0)
    # m.fs.permeate.split_fraction[0, "to_demin", "TDS"].fix(1.0)
    m.fs.permeate.split_fraction[0, "to_raw_water_tank_mix", "H2O"].fix(49 / 56.2)
    m.fs.permeate.split_fraction[0, "to_raw_water_tank_mix", "TDS"].fix(49 / 56.2)
    touch_flow_and_conc(m.fs.permeate)

    # ------------------------ COOLING TOWER ------------------------
    m.fs.cooling_tower = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["to_evaporation", "to_ww_surge_tank"],
        split_basis=SplittingType.componentFlow,
    )
    m.fs.cooling_tower.split_fraction[0, "to_evaporation", "H2O"].fix(9178.9 / 10522.4)
    m.fs.cooling_tower.split_fraction[0, "to_evaporation", "TDS"].fix(0)
    touch_flow_and_conc(m.fs.cooling_tower)

    m.fs.power_gen = Var(
        initialize=100, units=pyunits.MWh, bounds=(0, 200), doc="Power generated by SRP"
    )
    m.fs.power_gen_slope = Param(
        initialize=6.96333391e-03,
        mutable=True,
        units=pyunits.MWh / (pyunits.gallons / pyunits.minute),
        doc="Slope of power generated vs flow rate",
    )
    m.fs.power_gen_intercept = Param(
        initialize=28.89,
        mutable=True,
        units=pyunits.MWh,
        doc="Intercept of power generated vs flow rate",
    )
    m.fs.power_gen_constr = Constraint(
        expr=m.fs.power_gen
        == pyunits.convert(
            m.fs.power_gen_slope
            * pyunits.convert(
                m.fs.cooling_tower.to_evaporation_state[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallons / pyunits.minute,
            )
            + m.fs.power_gen_intercept,
            to_units=pyunits.MWh,
        )
    )

    # ------------------------ EVAPORATION FROM COOLING TOWER ------------------------
    m.fs.evaporation = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.evaporation)

    # ------------------------ WW SURGE TANK ------------------------
    m.fs.ww_surge_tank = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["to_ro_reject_tank", "to_bc_feed_tank", "to_mystery"],
        split_basis=SplittingType.componentFlow,
    )
    # original split from spreadsheet
    # modified so most of "mystery" missing flow goes to BC feed tank
    # so that the total flow to BCs is 350 + 350 + 450 gpm (the stated capacity of each BC)
    # m.fs.ww_surge_tank.split_fraction[0, "to_bc_feed_tank", "H2O"].fix(544.5 / 1353.1)
    # m.fs.ww_surge_tank.split_fraction[0, "to_bc_feed_tank", "TDS"].fix(544.5 / 1353.1)
    m.fs.ww_surge_tank.split_fraction[0, "to_bc_feed_tank", "H2O"].fix(966.2 / 1353.1)
    m.fs.ww_surge_tank.split_fraction[0, "to_bc_feed_tank", "TDS"].fix(966.2 / 1353.1)
    m.fs.ww_surge_tank.split_fraction[0, "to_ro_reject_tank", "H2O"].fix(305.1 / 1353.1)
    m.fs.ww_surge_tank.split_fraction[0, "to_ro_reject_tank", "TDS"].fix(305.1 / 1353.1)
    touch_flow_and_conc(m.fs.ww_surge_tank)

    # ------------------------ MYSTERY ------------------------
    m.fs.mystery = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.mystery)

    # ------------------------ BC FEED TANK ------------------------
    m.fs.bc_feed_tank = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["from_ro_reject_tank", "from_ww_surge_tank"],
    )
    m.fs.bc_feed_tank.outlet.pressure[0].fix(101325)
    touch_flow_and_conc(m.fs.bc_feed_tank)

    # ------------------------ RO REJECT TANK ------------------------
    m.fs.ro_reject_tank = Separator(
        property_package=m.fs.properties_feed,
        outlet_list=["to_bc_feed_tank", "to_ro_containment", "to_conc_waste", "to_uf"],
        split_basis=SplittingType.componentFlow,
    )
    m.fs.ro_reject_tank.split_fraction[0, "to_uf", "H2O"].fix(114.2 / 305.1)
    m.fs.ro_reject_tank.split_fraction[0, "to_uf", "TDS"].fix(114.2 / 305.1)
    m.fs.ro_reject_tank.split_fraction[0, "to_bc_feed_tank", "H2O"].fix(183.5 / 305.1)
    m.fs.ro_reject_tank.split_fraction[0, "to_bc_feed_tank", "TDS"].fix(183.5 / 305.1)
    m.fs.ro_reject_tank.split_fraction[0, "to_conc_waste", "H2O"].fix(6.4 / 305.1)
    m.fs.ro_reject_tank.split_fraction[0, "to_conc_waste", "TDS"].fix(6.4 / 305.1)
    touch_flow_and_conc(m.fs.ro_reject_tank)

    # ------------------------ RO CONTAINMENT ------------------------
    m.fs.ro_containment = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["from_uf", "from_ro_reject_tank", "from_ro"],
    )
    m.fs.ro_containment.outlet.pressure[0].fix(101325)
    touch_flow_and_conc(m.fs.ro_containment)

    # ------------------------ CONCENTRATED WASTE TANK ------------------------
    m.conc_waste_inlets = ["from_ro_reject_tank"]
    for bc_label in BCs:
        m.conc_waste_inlets.append(f"from_{bc_label}")
    m.fs.conc_waste = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=m.conc_waste_inlets,
    )
    touch_flow_and_conc(m.fs.conc_waste)

    # ------------------------ CONNECTIONS TO BCs ------------------------
    if len(BCs) > 1:
        # if there is more than one BC, use a separator to split flow to each
        m.BCs_outlets = []
        for bc_label in BCs:
            m.BCs_outlets.append(f"to_{bc_label}")
        m.fs.BCs = Separator(
            property_package=m.fs.properties_feed,
            outlet_list=m.BCs_outlets,
            split_basis=SplittingType.componentFlow,
        )
        # assume even split to each BC
        # even_split = 1 / len(BCs)
        split = 350 / (350 * 2 + 450)
        for i, bc_label in enumerate(BCs):
            if i == 0:
                continue
            m.fs.BCs.split_fraction[0, f"to_{bc_label}", "H2O"].fix(split)
            m.fs.BCs.split_fraction[0, f"to_{bc_label}", "TDS"].fix(split)
    else:
        m.fs.BCs = StateJunction(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.BCs)

    # ------------------------ EVAPORATION PONDS ------------------------
    m.fs.evaporation_ponds = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["from_conc_waste", "from_ro_containment"],
    )
    m.fs.evaporation_ponds.outlet.pressure[0].fix(101325)
    touch_flow_and_conc(m.fs.evaporation_ponds)

    # ------------------------ EVAPORATION PONDS ------------------------
    m.fs.pond_evap = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.pond_evap)

    m.demin_inlets = ["from_ro_permeate"]
    for bc_label in m.BCs:
        m.demin_inlets.append(f"from_{bc_label}")
    m.fs.demin = Mixer(
        property_package=m.fs.properties_feed,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=m.demin_inlets,
    )
    m.fs.demin.outlet.pressure[0].fix(101325)
    touch_flow_and_conc(m.fs.demin)

    m.fs.product = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.product)

    # =======================================================================
    # Connections
    # =======================================================================

    # WELLS TO RAW WATER TANK MIXER
    m.fs.wells_to_raw_water_tank_mix = Arc(
        source=m.fs.wells.outlet, destination=m.fs.raw_water_tank_mix.from_wells
    )
    m.fs.permeate_to_raw_water_tank_mix = Arc(
        source=m.fs.permeate.to_raw_water_tank_mix,
        destination=m.fs.raw_water_tank_mix.from_permeate,
    )

    # RAW WATER TANK MIXER TO RAW WATER TANK

    m.fs.raw_water_tank_mix_to_raw_water_tank = Arc(
        source=m.fs.raw_water_tank_mix.outlet, destination=m.fs.raw_water_tank.inlet
    )

    # RAW WATER TANK TO COOLING TOWER AND SERVICE & FIRE

    m.fs.raw_water_tank_to_cooling_tower = Arc(
        source=m.fs.raw_water_tank.to_cooling_tower,
        destination=m.fs.cooling_tower.inlet,
    )
    m.fs.raw_water_tank_to_service_and_fire = Arc(
        source=m.fs.raw_water_tank.to_service_and_fire,
        destination=m.fs.service_and_fire.inlet,
    )

    # COOLING TOWER TO EVAPORATION AND WW SURGE TANK
    m.fs.cooling_tower_to_evaporation = Arc(
        source=m.fs.cooling_tower.to_evaporation, destination=m.fs.evaporation.inlet
    )
    m.fs.cooling_tower_to_ww_surge_tank = Arc(
        source=m.fs.cooling_tower.to_ww_surge_tank, destination=m.fs.ww_surge_tank.inlet
    )

    # WW SURGE TANK TO RO REJECT TANK, BC FEED TANK, AND MYSTERY
    m.fs.ww_surge_tank_to_ro_reject_tank = Arc(
        source=m.fs.ww_surge_tank.to_ro_reject_tank,
        destination=m.fs.ro_reject_tank.inlet,
    )
    m.fs.ww_surge_tank_to_bc_feed_tank = Arc(
        source=m.fs.ww_surge_tank.to_bc_feed_tank,
        destination=m.fs.bc_feed_tank.from_ww_surge_tank,
    )
    m.fs.ww_surge_tank_to_mystery = Arc(
        source=m.fs.ww_surge_tank.to_mystery, destination=m.fs.mystery.inlet
    )

    # RO REJECT TANK TO BC FEED TANK, RO CONTAINMENT, CONC WASTE, AND UF
    m.fs.ro_reject_tank_to_bc_feed_tank = Arc(
        source=m.fs.ro_reject_tank.to_bc_feed_tank,
        destination=m.fs.bc_feed_tank.from_ro_reject_tank,
    )
    m.fs.ro_reject_tank_to_ro_containment = Arc(
        source=m.fs.ro_reject_tank.to_ro_containment,
        destination=m.fs.ro_containment.from_ro_reject_tank,
    )
    m.fs.ro_reject_tank_to_uf = Arc(
        source=m.fs.ro_reject_tank.to_uf, destination=m.fs.uf.inlet
    )
    m.fs.ro_reject_tank_to_conc_waste = Arc(
        source=m.fs.ro_reject_tank.to_conc_waste,
        destination=m.fs.conc_waste.from_ro_reject_tank,
        # destination=m.fs.conc_waste.inlet,
    )

    # UF TO RO CONTAINMENT AND RO PUMP

    m.fs.uf_to_ro_containment = Arc(
        source=m.fs.uf.to_ro_containment, destination=m.fs.ro_containment.from_uf
    )
    m.fs.uf_to_ro_pump = Arc(source=m.fs.uf.to_ro_pump, destination=m.fs.ro_pump.inlet)

    # RO PUMP TO RO
    m.fs.ro_pump_to_ro = Arc(source=m.fs.ro_pump.outlet, destination=m.fs.ro.inlet)

    # RO TO RO CONTAINMENT AND RO PERMEATE
    m.fs.ro_to_ro_containment = Arc(
        source=m.fs.ro.to_ro_containment, destination=m.fs.ro_containment.from_ro
    )
    m.fs.ro_to_ro_permeate = Arc(
        source=m.fs.ro.to_ro_permeate, destination=m.fs.permeate.inlet
    )

    m.fs.conc_waste_to_evaporation_ponds = Arc(
        source=m.fs.conc_waste.outlet,
        destination=m.fs.evaporation_ponds.from_conc_waste,
    )

    m.fs.ro_containment_to_evaporation_ponds = Arc(
        source=m.fs.ro_containment.outlet,
        destination=m.fs.evaporation_ponds.from_ro_containment,
    )

    m.fs.evaporation_ponds_to_pond_evap = Arc(
        source=m.fs.evaporation_ponds.outlet, destination=m.fs.pond_evap.inlet
    )

    m.fs.bc_feed_tank_to_bcs = Arc(
        source=m.fs.bc_feed_tank.outlet, destination=m.fs.BCs.inlet
    )

    m.fs.ro_permeate_to_demin = Arc(
        source=m.fs.permeate.to_demin, destination=m.fs.demin.from_ro_permeate
    )

    m.fs.demin_to_product = Arc(
        source=m.fs.demin.outlet, destination=m.fs.product.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # =======================================================================
    # Scaling
    # =======================================================================

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

    iscale.calculate_scaling_factors(m)

    m.fs.wells.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin,
            ("conc_mass_phase_comp", ("Liq", "TDS")): Cin,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + feed_temp,
        },
        hold_state=True,
    )

    # =======================================================================
    # Initialization
    # =======================================================================

    print(f"dof = {degrees_of_freedom(m)}")
    m.fs.wells.initialize()
    propagate_state(m.fs.wells_to_raw_water_tank_mix)
    perm_flow = pyunits.convert(
        49 * pyunits.gallon / pyunits.minute, to_units=pyunits.m**3 / pyunits.s
    )()
    m.fs.raw_water_tank_mix.from_permeate_state.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): perm_flow,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 0.5,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + 27,
        },
        hold_state=False,
    )
    m.fs.raw_water_tank_mix.initialize()
    propagate_state(m.fs.raw_water_tank_mix_to_raw_water_tank)

    m.fs.raw_water_tank.initialize()
    propagate_state(m.fs.raw_water_tank_to_cooling_tower)
    propagate_state(m.fs.raw_water_tank_to_service_and_fire)

    m.fs.cooling_tower.initialize()
    propagate_state(m.fs.cooling_tower_to_evaporation)
    propagate_state(m.fs.cooling_tower_to_ww_surge_tank)
    m.fs.evaporation.initialize()

    m.fs.service_and_fire.initialize()

    m.fs.ro_pump.initialize()
    propagate_state(m.fs.ro_pump_to_ro)

    m.fs.ro.initialize()
    propagate_state(m.fs.ro_to_ro_containment)
    propagate_state(m.fs.ro_to_ro_permeate)

    m.fs.permeate.initialize()
    propagate_state(m.fs.permeate_to_raw_water_tank_mix)

    m.fs.ww_surge_tank.initialize()
    propagate_state(m.fs.ww_surge_tank_to_ro_reject_tank)
    propagate_state(m.fs.ww_surge_tank_to_bc_feed_tank)
    propagate_state(m.fs.ww_surge_tank_to_mystery)

    m.fs.mystery.initialize()

    m.fs.ro_reject_tank.initialize()
    propagate_state(m.fs.ro_reject_tank_to_ro_containment)
    propagate_state(m.fs.ro_reject_tank_to_bc_feed_tank)
    propagate_state(m.fs.ro_reject_tank_to_conc_waste)
    propagate_state(m.fs.ro_reject_tank_to_uf)

    m.fs.uf.initialize()
    propagate_state(m.fs.uf_to_ro_containment)
    propagate_state(m.fs.uf_to_ro_pump)

    m.fs.bc_feed_tank.initialize()
    propagate_state(m.fs.bc_feed_tank_to_bcs)

    m.fs.BCs.initialize()

    m.fs.conc_waste.initialize()
    propagate_state(m.fs.conc_waste_to_evaporation_ponds)

    m.fs.ro_containment.initialize()
    propagate_state(m.fs.conc_waste_to_evaporation_ponds)
    propagate_state(m.fs.ro_containment_to_evaporation_ponds)

    propagate_state(m.fs.ro_permeate_to_demin)
    m.fs.demin.initialize()
    propagate_state(m.fs.demin_to_product)

    m.fs.product.initialize()

    m.fs.evaporation_ponds.initialize()
    propagate_state(m.fs.evaporation_ponds_to_pond_evap)

    m.fs.pond_evap.initialize()

    print(f"dof = {degrees_of_freedom(m)}")

    return m


def add_MVCs(m, recovery_vol=0.92):

    for bc_label in m.BCs:
        bc_fs = FlowsheetBlock()
        m.fs.add_component(bc_label, bc_fs)
        bc_fs = m.fs.find_component(bc_label)
        bc.build_mvc(m, bc_fs)
        bc.add_mvc_costing(m, bc_fs)

    m.fs.costing.cost_process()

    for bc_label in m.BCs:
        bc_fs = m.fs.find_component(bc_label)
        print(f"dof {bc_label} = {degrees_of_freedom(bc_fs)}")
        bc.set_mvc_operating_conditions(m, bc_fs)
        bc.set_mvc_scaling(m, bc_fs)
        if len(m.BCs) > 1:
            feed_props = m.fs.BCs.find_component(f"to_{bc_label}_state")
            upstream = m.fs.BCs.find_component(f"to_{bc_label}")
        else:
            feed_props = m.fs.BCs.find_component("properties")
            upstream = m.fs.BCs.find_component("outlet")
        upstream_to_mvc = Arc(source=upstream, destination=bc_fs.feed.inlet)
        bc.init_mvc(m, bc_fs, feed_props=feed_props[0])
        bc_fs.add_component(f"upstream_to_{bc_label}", upstream_to_mvc)
        propagate_state(upstream_to_mvc)
        TransformationFactory("network.expand_arcs").apply_to(m)
        bc_fs.feed.initialize()
        bc_fs.recovery_mass.unfix()
        bc_fs.recovery_vol.fix(recovery_vol)
        results = bc.solve_mvc(bc_fs)
        bc_fs.feed.properties[0].flow_vol_phase
        bc_fs.feed.properties[0].conc_mass_phase_comp
        bc_fs.feed.initialize()
        bc_fs.disposal.properties[0].flow_vol_phase
        bc_fs.disposal.properties[0].conc_mass_phase_comp
        bc_fs.disposal.initialize()
        bc_fs.product.properties[0].flow_vol_phase
        bc_fs.product.properties[0].conc_mass_phase_comp
        bc_fs.product.initialize()
        # results = bc.solve_m÷vc(bc_fs)

    results = bc.solve_mvc(m)

    for bc_label in m.BCs:
        bc_fs = m.fs.find_component(bc_label)
        print(f"dof {bc_label} = {degrees_of_freedom(bc_fs)}")
        bc_fs.compressor.pressure_ratio.fix(1.6)

    # m.fs.obj = Objective(expr=m.fs.costing.SEC)
    print(f"dof = {degrees_of_freedom(m)}")

    results = bc.solve_mvc(m)

    for bc_label in m.BCs:
        bc_fs = m.fs.find_component(bc_label)
        conc_waste_port = m.fs.conc_waste.find_component(f"from_{bc_label}")
        bc_to_conc_waste = Arc(
            source=bc_fs.disposal.outlet, destination=conc_waste_port
        )
        m.fs.add_component(f"{bc_label}_to_conc_waste", bc_to_conc_waste)
        bc_to_conc_waste = m.fs.find_component(f"{bc_label}_to_conc_waste")
        propagate_state(bc_to_conc_waste)

    for bc_label in m.BCs:
        bc_fs = m.fs.find_component(bc_label)
        demin_port = m.fs.demin.find_component(f"from_{bc_label}")
        bc_to_demin = Arc(source=bc_fs.product.outlet, destination=demin_port)
        m.fs.add_component(f"{bc_label}_to_demin", bc_to_demin)
        bc_to_demin = m.fs.find_component(f"{bc_label}_to_demin")
        propagate_state(bc_to_demin)

    TransformationFactory("network.expand_arcs").apply_to(m)
    print(f"dof = {degrees_of_freedom(m)}")

    iscale.calculate_scaling_factors(m)

    d = {"evaporator.area": 455.0, "hx_brine.area": 10.4, "hx_distillate.area": 60.8}

    for i, bc_label in enumerate(m.BCs):
        if i == 0:
            bc0 = m.fs.find_component(bc_label)
            # for k, v in d.items():
            #     a = bc0.find_component(k)
            #     # print(a.name)
            #     # if a is not None:
            #     a.setlb(0.95 * v)
            #     a.setub(1.05 * v)
            #     a.set_value(v)
            continue
        bc_fs = m.fs.find_component(bc_label)
        evap_area_constr = Constraint(expr=bc0.evaporator.area == bc_fs.evaporator.area)
        hx_brine_area_constr = Constraint(expr=bc0.hx_brine.area == bc_fs.hx_brine.area)
        hx_distillate_area_constr = Constraint(
            expr=bc0.hx_distillate.area == bc_fs.hx_distillate.area
        )
        # for k, v in d.items():
        #     a = bc_fs.find_component(k)
        #     # print(a.name)
        #     # if a is not None:
        #     a.setlb(0.95 * v)
        #     a.setub(1.05 * v)
        #     a.set_value(v)
        m.fs.add_component(f"{bc_label}_evap_area_constr", evap_area_constr)
        m.fs.add_component(f"{bc_label}_hx_brine_area_constr", hx_brine_area_constr)
        m.fs.add_component(
            f"{bc_label}_hx_distillate_area_constr", hx_distillate_area_constr
        )
    # assert False
    print(f"dof = {degrees_of_freedom(m)}")
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
    )


def print_stream_flows(m, w=30):

    for b in m.fs.component_objects(Block, descend_into=False):
        if "_expanded" in b.name:
            continue
        if b is m.fs.properties_feed:
            continue
        if b is m.fs.properties_vapor:
            continue

        # title = "KBHDP RPT 1 Report"
        title = b.name.replace("fs.", "").replace("_", " ").upper()
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "=" * side + f" {title} " + "=" * side
        print(f"\n{header}\n")
        if (
            isinstance(b, Feed)
            or isinstance(b, Product)
            or isinstance(b, StateJunction)
        ):
            props = b.find_component("properties")
            flow_in = value(
                pyunits.convert(
                    props[0].flow_vol_phase["Liq"],
                    to_units=pyunits.gallons / pyunits.minute,
                )
            )
            conc_in = value(
                pyunits.convert(
                    props[0].conc_mass_phase_comp["Liq", "TDS"],
                    to_units=pyunits.mg / pyunits.L,
                )
            )
            print(f'{"Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
            print(f'{"TDS":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
        elif isinstance(b, Separator) or isinstance(b, Mixer):
            if b.name == "fs.ro":
                recov = b.split_fraction[0, "to_ro_permeate", "H2O"]() * 100
                print(f'{"Recovery":<{w}s}{f"{recov:<{w},.2f}"}{"%":<{w}s}')
            ms = b.find_component("mixed_state")
            if isinstance(b, Separator):
                flow_in = value(
                    pyunits.convert(
                        ms[0].flow_vol_phase["Liq"],
                        to_units=pyunits.gallons / pyunits.minute,
                    )
                )
                conc_in = value(
                    pyunits.convert(
                        ms[0].conc_mass_phase_comp["Liq", "TDS"],
                        to_units=pyunits.mg / pyunits.L,
                    )
                )
                print(f'{"Feed Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
                print(f'{"Feed TDS":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
                tot_flow_out = sum(
                    value(
                        value(
                            pyunits.convert(
                                b.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                                to_units=pyunits.gallons / pyunits.minute,
                            )
                        )
                    )
                    for x in b.config.outlet_list
                )
                print(
                    f'{"TOTAL OUTLET FLOW":<{w}s}{f"{tot_flow_out:<{w},.1f}"}{"gpm":<{w}s}'
                )
                for x in b.config.outlet_list:
                    sb = b.find_component(f"{x}_state")
                    flow_out = value(
                        pyunits.convert(
                            sb[0].flow_vol_phase["Liq"],
                            to_units=pyunits.gallons / pyunits.minute,
                        )
                    )
                    conc_out = value(
                        pyunits.convert(
                            sb[0].conc_mass_phase_comp["Liq", "TDS"],
                            to_units=pyunits.mg / pyunits.L,
                        )
                    )
                    print(
                        f'{"   Flow " + x.replace("_", " ").title():<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}'
                    )
                    print(
                        f'{"   TDS " + x.replace("_", " ").title():<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}'
                    )
            elif isinstance(b, Mixer):
                tot_flow_in = sum(
                    value(
                        value(
                            pyunits.convert(
                                b.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                                to_units=pyunits.gallons / pyunits.minute,
                            )
                        )
                    )
                    for x in b.config.inlet_list
                )
                print(
                    f'{"TOTAL INLET FLOW":<{w}s}{f"{tot_flow_in:<{w},.1f}"}{"gpm":<{w}s}'
                )
                for x in b.config.inlet_list:
                    sb = b.find_component(f"{x}_state")
                    flow_in = value(
                        pyunits.convert(
                            sb[0].flow_vol_phase["Liq"],
                            to_units=pyunits.gallons / pyunits.minute,
                        )
                    )
                    tot_flow_in += flow_in
                    conc_in = value(
                        pyunits.convert(
                            sb[0].conc_mass_phase_comp["Liq", "TDS"],
                            to_units=pyunits.mg / pyunits.L,
                        )
                    )
                    print(
                        f'{"   Flow " + x.replace("_", " ").title():<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}'
                    )
                    print(
                        f'{"   TDS " + x.replace("_", " ").title():<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}'
                    )
                flow_out = value(
                    pyunits.convert(
                        ms[0].flow_vol_phase["Liq"],
                        to_units=pyunits.gallons / pyunits.minute,
                    )
                )
                conc_out = value(
                    pyunits.convert(
                        ms[0].conc_mass_phase_comp["Liq", "TDS"],
                        to_units=pyunits.mg / pyunits.L,
                    )
                )
                print(f'{"Outlet Flow":<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}')
                print(f'{"Outlet TDS":<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}')
        elif isinstance(b, Pump):
            cv = b.find_component("control_volume")
            flow_in = value(
                pyunits.convert(
                    cv.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.gallons / pyunits.minute,
                )
            )
            conc_in = value(
                pyunits.convert(
                    cv.properties_in[0].conc_mass_phase_comp["Liq", "TDS"],
                    to_units=pyunits.mg / pyunits.L,
                )
            )
            pressure = value(
                pyunits.convert(cv.properties_in[0].pressure, to_units=pyunits.bar)
            )
            print(f'{"Inlet Pressure":<{w}s}{f"{pressure:<{w},.1f}"}{"bar":<{w}s}')
            print(f'{"Inlet Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
            print(f'{"Inlet TDS":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
            flow_out = value(
                pyunits.convert(
                    cv.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.gallons / pyunits.minute,
                )
            )
            conc_out = value(
                pyunits.convert(
                    cv.properties_out[0].conc_mass_phase_comp["Liq", "TDS"],
                    to_units=pyunits.mg / pyunits.L,
                )
            )
            pressure = value(
                pyunits.convert(cv.properties_out[0].pressure, to_units=pyunits.bar)
            )
            print(f'{"Outlet Flow":<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}')
            print(f'{"Outlet TDS":<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}')
            print(f'{"Outlet Pressure":<{w}s}{f"{pressure:<{w},.1f}"}{"bar":<{w}s}')
        elif isinstance(b, FlowsheetBlock):
            print_MVC_stream_flows(b, w=w)


def touch_flow_and_conc(b):
    props = b.find_component("properties")
    if props is not None:
        props[0].flow_vol_phase
        props[0].conc_mass_phase_comp
        return
    cv = b.find_component("control_volume")
    if cv is not None:
        cv.properties_in[0].flow_vol_phase
        cv.properties_in[0].conc_mass_phase_comp
        cv.properties_out[0].flow_vol_phase
        cv.properties_out[0].conc_mass_phase_comp
        return
    b.mixed_state[0].flow_vol_phase
    b.mixed_state[0].conc_mass_phase_comp
    if isinstance(b, Mixer):
        for x in b.config.inlet_list:
            sb = b.find_component(f"{x}_state")
            sb[0].flow_vol_phase
            sb[0].conc_mass_phase_comp
    if isinstance(b, Separator):
        for x in b.config.outlet_list:
            sb = b.find_component(f"{x}_state")
            sb[0].flow_vol_phase
            sb[0].conc_mass_phase_comp


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

    m = build_primary_fs()
    m.fs.ro_pump.control_volume.properties_out[0].pressure.unfix()
    m.fs.ro.split_fraction[0, "to_ro_permeate", "H2O"].fix(0.7)
    results = solver.solve(m)
    print_infeasible_constraints(m)
    print(f"termination condition: {results.solver.termination_condition}")
    # bc_fs.add_MVCs(m)
    print_stream_flows(m)
    add_MVCs(m, recovery_vol=0.92)
