from pyomo.environ import *
from pyomo.environ import units as pyunits
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Separator, Mixer, Product, StateJunction
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.separator import SplittingType
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap.costing import WaterTAPCosting
from watertap.core.solvers import get_solver

from srp.components import *
from srp.utils.utils import touch_flow_and_conc


def build_srp(
    Qin=11343, Cin=1467, feed_temp=27, BCs=["BC_A", "BC_B", "BC_C"], perm_flow_guess=49
):

    Qin = Qin * pyunits.gallons / pyunits.minute
    Cin = Cin * pyunits.mg / pyunits.L
    perm_flow_guess = perm_flow_guess * pyunits.gallons / pyunits.minute

    m = ConcreteModel()
    m.Qin = Qin
    m.Cin = Cin
    m.feed_temp = feed_temp
    m.perm_flow_guess = perm_flow_guess
    m.BCs = BCs
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.costing = WaterTAPCosting()

    m.fs.properties_feed = SeawaterParameterBlock()
    m.fs.properties_vapor = SteamParameterBlock()

    # ------------------------ WELLS ------------------------
    m.fs.wells = Feed(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.wells)

    # ------------------------ RAW WATER TANK MIXER ------------------------

    m.fs.raw_water_tank_mixer = FlowsheetBlock(dynamic=False)
    build_mixer(
        m.fs.raw_water_tank_mixer,
        prop_package=m.fs.properties_feed,
        inlet_list=["from_permeate", "from_wells"],
    )

    # ------------------------ RAW WATER TANK ------------------------

    m.fs.raw_water_tank = FlowsheetBlock(dynamic=False)
    build_separator(
        m.fs.raw_water_tank,
        outlet_list=["to_cooling_tower", "to_service_and_fire"],
        prop_package=m.fs.properties_feed,
    )

    # ------------------------ SERVICE AND FIRE ------------------------
    m.fs.service_and_fire = FlowsheetBlock(dynamic=False)
    build_product(m.fs.service_and_fire, prop_package=m.fs.properties_feed)

    # ------------------------ UF SYSTEM ------------------------
    m.fs.uf = FlowsheetBlock(dynamic=False)
    outlet_list = ["to_ro_pump", "to_ro_containment"]
    build_separator(m.fs.uf, outlet_list=outlet_list, prop_package=m.fs.properties_feed)

    # ------------------------ RO PUMP ------------------------

    m.fs.ro_pump = FlowsheetBlock(dynamic=False)
    build_pump(m.fs.ro_pump, prop_package=m.fs.properties_feed)

    # ------------------------ RO ------------------------

    m.fs.ro = FlowsheetBlock(dynamic=False)
    build_ro(
        m.fs.ro,
        pump=m.fs.ro_pump.unit,
        outlet_list=["to_ro_containment", "to_ro_permeate"],
        prop_package=m.fs.properties_feed,
    )

    # ------------------------ PERMEATE ------------------------

    m.fs.permeate = FlowsheetBlock(dynamic=False)
    build_separator(
        m.fs.permeate,
        prop_package=m.fs.properties_feed,
        outlet_list=["to_demin", "to_raw_water_tank_mix"],
    )

    # ------------------------ COOLING TOWER ------------------------

    m.fs.cooling_tower = FlowsheetBlock(dynamic=False)
    build_separator(
        m.fs.cooling_tower,
        prop_package=m.fs.properties_feed,
        outlet_list=["to_evaporation", "to_ww_surge_tank"],
    )

    # ------------------------ EVAPORATION FROM COOLING TOWER ------------------------
    m.fs.cooling_tower_evap = FlowsheetBlock(dynamic=False)
    build_product(m.fs.cooling_tower_evap, prop_package=m.fs.properties_feed)

    # ------------------------ WW SURGE TANK ------------------------

    m.fs.ww_surge_tank = FlowsheetBlock(dynamic=False)
    build_separator(
        m.fs.ww_surge_tank,
        prop_package=m.fs.properties_feed,
        outlet_list=["to_ro_reject_tank", "to_bc_feed_tank", "to_mystery"],
    )

    # ------------------------ MYSTERY ------------------------
    m.fs.mystery = FlowsheetBlock(dynamic=False)
    build_product(m.fs.mystery, prop_package=m.fs.properties_feed)

    # ------------------------ BC FEED TANK ------------------------

    m.fs.bc_feed_tank = FlowsheetBlock(dynamic=False)
    build_mixer(
        m.fs.bc_feed_tank,
        prop_package=m.fs.properties_feed,
        inlet_list=["from_ro_reject_tank", "from_ww_surge_tank"],
    )
    # ------------------------ RO REJECT TANK ------------------------

    m.fs.ro_reject_tank = FlowsheetBlock(dynamic=False)
    build_separator(
        m.fs.ro_reject_tank,
        prop_package=m.fs.properties_feed,
        outlet_list=["to_bc_feed_tank", "to_ro_containment", "to_conc_waste", "to_uf"],
    )
    # ------------------------ RO CONTAINMENT ------------------------

    m.fs.ro_containment = FlowsheetBlock(dynamic=False)
    build_mixer(
        m.fs.ro_containment,
        prop_package=m.fs.properties_feed,
        inlet_list=["from_uf", "from_ro_reject_tank", "from_ro"],
    )

    # ------------------------ CONCENTRATED WASTE TANK ------------------------
    m.fs.conc_waste = FlowsheetBlock(dynamic=False)
    build_mixer(
        m.fs.conc_waste,
        prop_package=m.fs.properties_feed,
        inlet_list=["from_ro_reject_tank"],
    )

    # ------------------------ CONNECTIONS TO BCs ------------------------
    if len(m.BCs) > 0:
        m.BC_outlets = []
        for bc_label in m.BCs:
            m.BC_outlets.append(f"to_{bc_label.lower()}")
        print(m.BC_outlets)
        # assert False
        m.fs.BCs = FlowsheetBlock(dynamic=False)
        build_separator(
            m.fs.BCs, prop_package=m.fs.properties_feed, outlet_list=m.BC_outlets
        )
        split = 350 / (350 * 2 + 450)
        for i, bc_label in enumerate(m.BC_outlets):
            if i == 0:
                continue
            m.fs.BCs.unit.split_fraction[0, f"{bc_label}", "H2O"].fix(split)
            m.fs.BCs.unit.split_fraction[0, f"{bc_label}", "TDS"].fix(split)

    # ------------------------ EVAPORATION PONDS ------------------------
    m.fs.evaporation_ponds = FlowsheetBlock(dynamic=False)
    build_mixer(
        m.fs.evaporation_ponds,
        prop_package=m.fs.properties_feed,
        inlet_list=["from_conc_waste", "from_ro_containment"],
    )

    # ------------------------ POND EVAPORATION ------------------------
    m.fs.pond_evap = FlowsheetBlock(dynamic=False)
    build_product(m.fs.pond_evap, prop_package=m.fs.properties_feed)

    # ------------------------ DEMIN ------------------------
    m.demin_inlets = ["from_ro_permeate"]
    for bc_label in m.BCs:
        m.demin_inlets.append(f"from_{bc_label.lower()}")
    m.fs.demin = FlowsheetBlock(dynamic=False)
    build_mixer(
        m.fs.demin, prop_package=m.fs.properties_feed, inlet_list=m.demin_inlets
    )

    m.fs.product = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.product)

    return m


def connect_srp(m):

    # WELLS TO RAW WATER TANK MIXER
    m.fs.wells_to_raw_water_tank_mix = Arc(
        source=m.fs.wells.outlet, destination=m.fs.raw_water_tank_mixer.unit.from_wells
    )
    # PERMEATE TO RAW WATER TANK MIXER
    m.fs.permeate_to_raw_water_tank_mix = Arc(
        source=m.fs.permeate.unit.to_raw_water_tank_mix,
        destination=m.fs.raw_water_tank_mixer.unit.from_permeate,
    )

    # RAW WATER TANK MIXER TO RAW WATER TANK

    m.fs.raw_water_tank_mix_to_raw_water_tank = Arc(
        source=m.fs.raw_water_tank_mixer.product.outlet,
        destination=m.fs.raw_water_tank.feed.inlet,
    )

    # RAW WATER TANK TO COOLING TOWER AND SERVICE & FIRE

    m.fs.raw_water_tank_to_cooling_tower = Arc(
        source=m.fs.raw_water_tank.to_cooling_tower.outlet,
        destination=m.fs.cooling_tower.feed.inlet,
    )
    m.fs.raw_water_tank_to_service_and_fire = Arc(
        source=m.fs.raw_water_tank.to_service_and_fire.outlet,
        destination=m.fs.service_and_fire.feed.inlet,
    )

    # COOLING TOWER TO EVAPORATION AND WW SURGE TANK
    m.fs.cooling_tower_to_evaporation = Arc(
        source=m.fs.cooling_tower.to_evaporation.outlet,
        destination=m.fs.cooling_tower_evap.feed.inlet,
    )
    m.fs.cooling_tower_to_ww_surge_tank = Arc(
        source=m.fs.cooling_tower.to_ww_surge_tank.outlet,
        destination=m.fs.ww_surge_tank.feed.inlet,
    )

    # WW SURGE TANK TO RO REJECT TANK, BC FEED TANK, AND MYSTERY
    m.fs.ww_surge_tank_to_ro_reject_tank = Arc(
        source=m.fs.ww_surge_tank.to_ro_reject_tank.outlet,
        destination=m.fs.ro_reject_tank.feed.inlet,
    )
    m.fs.ww_surge_tank_to_bc_feed_tank = Arc(
        source=m.fs.ww_surge_tank.to_bc_feed_tank.outlet,
        destination=m.fs.bc_feed_tank.from_ww_surge_tank.inlet,
    )
    m.fs.ww_surge_tank_to_mystery = Arc(
        source=m.fs.ww_surge_tank.to_mystery.outlet, destination=m.fs.mystery.feed.inlet
    )

    # RO REJECT TANK TO BC FEED TANK, RO CONTAINMENT, CONC WASTE, AND UF
    m.fs.ro_reject_tank_to_bc_feed_tank = Arc(
        source=m.fs.ro_reject_tank.to_bc_feed_tank.outlet,
        destination=m.fs.bc_feed_tank.from_ro_reject_tank.inlet,
    )
    m.fs.ro_reject_tank_to_ro_containment = Arc(
        source=m.fs.ro_reject_tank.to_ro_containment.outlet,
        destination=m.fs.ro_containment.from_ro_reject_tank.inlet,
    )
    m.fs.ro_reject_tank_to_uf = Arc(
        source=m.fs.ro_reject_tank.to_uf.outlet, destination=m.fs.uf.feed.inlet
    )
    m.fs.ro_reject_tank_to_conc_waste = Arc(
        source=m.fs.ro_reject_tank.to_conc_waste.outlet,
        destination=m.fs.conc_waste.from_ro_reject_tank.inlet,
        # destination=m.fs.conc_waste.inlet,
    )

    # UF TO RO CONTAINMENT AND RO PUMP

    m.fs.uf_to_ro_containment = Arc(
        source=m.fs.uf.to_ro_containment.outlet,
        destination=m.fs.ro_containment.from_uf.inlet,
    )
    m.fs.uf_to_ro_pump = Arc(
        source=m.fs.uf.to_ro_pump.outlet, destination=m.fs.ro_pump.feed.inlet
    )

    # RO PUMP TO RO
    m.fs.ro_pump_to_ro = Arc(
        source=m.fs.ro_pump.product.outlet, destination=m.fs.ro.feed.inlet
    )

    # RO TO RO CONTAINMENT AND RO PERMEATE
    m.fs.ro_to_ro_containment = Arc(
        source=m.fs.ro.to_ro_containment.outlet,
        destination=m.fs.ro_containment.from_ro.inlet,
    )
    m.fs.ro_to_ro_permeate = Arc(
        source=m.fs.ro.to_ro_permeate.outlet, destination=m.fs.permeate.feed.inlet
    )

    m.fs.conc_waste_to_evaporation_ponds = Arc(
        source=m.fs.conc_waste.product.outlet,
        destination=m.fs.evaporation_ponds.from_conc_waste.inlet,
    )

    m.fs.ro_containment_to_evaporation_ponds = Arc(
        source=m.fs.ro_containment.product.outlet,
        destination=m.fs.evaporation_ponds.from_ro_containment.inlet,
    )

    m.fs.evaporation_ponds_to_pond_evap = Arc(
        source=m.fs.evaporation_ponds.product.outlet,
        destination=m.fs.pond_evap.feed.inlet,
    )

    m.fs.bc_feed_tank_to_bcs = Arc(
        source=m.fs.bc_feed_tank.product.outlet, destination=m.fs.BCs.feed.inlet
    )

    m.fs.ro_permeate_to_demin = Arc(
        source=m.fs.permeate.to_demin.outlet,
        destination=m.fs.demin.from_ro_permeate.inlet,
    )

    m.fs.demin_to_product = Arc(
        source=m.fs.demin.product.outlet, destination=m.fs.product.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_srp_scaling(m):

    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s)() * 1000),
        index=("Liq", "H2O"),
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp",
        1 / (pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s)()),
        index=("Liq", "TDS"),
    )
    iscale.calculate_scaling_factors(m)


def set_srp_operating_conditions(m):

    m.fs.wells.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): m.Qin,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.Cin,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + m.feed_temp,
        },
        hold_state=True,
    )

    raw_water_tank_splits = {
        "to_cooling_tower": {"H2O": 10522.4 / 11343, "TDS": 10522.4 / 11343}
    }
    set_separator_op_conditions(
        m.fs.raw_water_tank, split_fractions=raw_water_tank_splits
    )

    uf_splits = {
        "to_ro_pump": {"H2O": 91.3 / 114.2, "TDS": 91.3 / 114.2}
    }
    set_separator_op_conditions(
        m.fs.uf, split_fractions=uf_splits
    )

    permeate_splits = {
        "to_raw_water_tank_mix": {"H2O": 49 / 56.2, "TDS": 49 / 56.2}
    }
    set_separator_op_conditions(
        m.fs.permeate, split_fractions=permeate_splits
    )

    cooling_tower_splits = {
        "to_evaporation": {"H2O": 9178.9 / 10522.4, "TDS": 0}
    }
    set_separator_op_conditions(
        m.fs.cooling_tower, split_fractions=cooling_tower_splits
    )

    ww_surge_tank_splits = {
        "to_bc_feed_tank": {"H2O": 966.2 / 1353.1, "TDS": 966.2 / 1353.1},
        "to_ro_reject_tank": {"H2O": 305.1 / 1353.1, "TDS": 305.1 / 1353.1},
    }
    set_separator_op_conditions(
        m.fs.ww_surge_tank, split_fractions=ww_surge_tank_splits
    )

    ro_rej_tank_splits = {
        "to_uf": {"H2O": 114.2 / 305.1, "TDS": 114.2 / 305.1},
        "to_bc_feed_tank": {"H2O": 183.5 / 305.1, "TDS": 183.5 / 305.1},
        "to_conc_waste": {"H2O": 6.4 / 305.1, "TDS": 6.4 / 305.1},
    }
    set_separator_op_conditions(
        m.fs.ro_reject_tank, split_fractions=ro_rej_tank_splits
    )

def initialize_srp(m):

    m.fs.wells.initialize()
    propagate_state(m.fs.wells_to_raw_water_tank_mix)

    m.fs.raw_water_tank_mixer.from_permeate.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): m.perm_flow_guess,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 0.5,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + 27,
        },
        hold_state=False,
    )

    m.fs.raw_water_tank_mixer.unit.from_permeate_state.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): m.perm_flow_guess,
            ("conc_mass_phase_comp", ("Liq", "TDS")): 0.5,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + 27,
        },
        hold_state=False,
    )

    init_mixer(m.fs.raw_water_tank_mixer)

    propagate_state(m.fs.raw_water_tank_mix_to_raw_water_tank)

    init_separator(m.fs.raw_water_tank)
    propagate_state(m.fs.raw_water_tank_to_cooling_tower)

    propagate_state(m.fs.raw_water_tank_to_service_and_fire)
    init_product(m.fs.service_and_fire)

    init_separator(m.fs.cooling_tower)
    propagate_state(m.fs.cooling_tower_to_evaporation)

    init_product(m.fs.cooling_tower_evap)
    propagate_state(m.fs.cooling_tower_to_ww_surge_tank)


def run_srp():

    m = build_srp()
    connect_srp(m)
    set_srp_scaling(m)
    set_srp_operating_conditions(m)
    initialize_srp(m)

    # assert degrees_of_freedom(m) == 0
    print(f"\nDegrees of freedom: {degrees_of_freedom(m)}\n")

    solver = get_solver()
    # results = solver.solve(m, tee=True)
    # assert_optimal_termination(results)


if __name__ == "__main__":
    run_srp()
