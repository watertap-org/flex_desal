from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Objective,
    Block,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Separator, Mixer, Product, StateJunction
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap.costing import WaterTAPCosting
from watertap.core.solvers import get_solver
from watertap.unit_models import Pump

from srp.components import *
from srp.utils.utils import touch_flow_and_conc

__all__ = [
    "build_srp",
    "connect_srp",
    "set_srp_scaling",
    "set_srp_operating_conditions",
    "add_bcs",
    "initialize_srp",
    "print_stream_flows",
    "add_bcs_basic",
    "init_bcs_basic",
    "add_bcs",
]

solver = get_solver()


def build_srp(
    Qin=11343,
    Cin=1467,
    feed_temp=27,
    BCs=["BC_A", "BC_B", "BC_C"],
    perm_flow_guess=49,
    add_basic_bcs=False,
):
    # BCs = []
    Qin = Qin * pyunits.gallons / pyunits.minute
    Cin = Cin * pyunits.mg / pyunits.L
    perm_flow_guess = perm_flow_guess * pyunits.gallons / pyunits.minute

    m = ConcreteModel()
    m.Qin = Qin
    m.Cin = Cin
    m.feed_temp = feed_temp
    m.perm_flow_guess = perm_flow_guess
    m.BCs = BCs
    m.add_basic_bcs = add_basic_bcs

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
    m.conc_waste_inlets = ["from_ro_reject_tank"]
    if len(m.BCs) > 0:
        for bc_label in m.BCs:
            m.conc_waste_inlets.append(f"from_{bc_label.lower()}")
    else:
        pass

    build_mixer(
        m.fs.conc_waste,
        prop_package=m.fs.properties_feed,
        inlet_list=m.conc_waste_inlets,
    )

    # ------------------------ CONNECTIONS TO BCs ------------------------
    m.fs.bcs = FlowsheetBlock(dynamic=False)
    if len(m.BCs) > 0:
        m.bc_outlets = []
        for bc_label in m.BCs:
            m.bc_outlets.append(f"to_{bc_label.lower()}")

        build_separator(
            m.fs.bcs, prop_package=m.fs.properties_feed, outlet_list=m.bc_outlets
        )
    else:
        m.fs.bcs.unit = StateJunction(property_package=m.fs.properties_feed)

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

    # ------------------------ FINAL PRODUCT ------------------------
    m.fs.product = Product(property_package=m.fs.properties_feed)
    touch_flow_and_conc(m.fs.product)

    if add_basic_bcs:
        add_bcs_basic(m)
    # else:
    #     add_bcs(m)

    # m.fs.costing.cost_process()

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
    if len(m.BCs) > 0:

        m.fs.bc_feed_tank_to_bcs = Arc(
            source=m.fs.bc_feed_tank.product.outlet, destination=m.fs.bcs.feed.inlet
        )
    else:

        m.fs.bc_feed_tank_to_bcs = Arc(
            source=m.fs.bc_feed_tank.product.outlet, destination=m.fs.bcs.unit.inlet
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

    set_pump_scaling(m.fs.ro_pump)

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

    uf_splits = {"to_ro_pump": {"H2O": 91.3 / 114.2, "TDS": 91.3 / 114.2}}
    set_separator_op_conditions(m.fs.uf, split_fractions=uf_splits)

    permeate_splits = {"to_raw_water_tank_mix": {"H2O": 49 / 56.2, "TDS": 49 / 56.2}}
    set_separator_op_conditions(m.fs.permeate, split_fractions=permeate_splits)

    cooling_tower_splits = {"to_evaporation": {"H2O": 9178.9 / 10522.4, "TDS": 0}}
    set_separator_op_conditions(
        m.fs.cooling_tower, split_fractions=cooling_tower_splits
    )

    ww_surge_tank_splits = {
        "to_bc_feed_tank": {"H2O": 966.2 / 1353.1, "TDS": 966.2 / 1353.1},
        "to_ro_reject_tank": {"H2O": 305.1 / 1353.1, "TDS": 305.1 / 1353.1},
        # "to_mystery": {"H2O": 0, "TDS": 0},
    }
    set_separator_op_conditions(
        m.fs.ww_surge_tank, split_fractions=ww_surge_tank_splits
    )

    ro_rej_tank_splits = {
        "to_uf": {"H2O": 114.2 / 305.1, "TDS": 114.2 / 305.1},
        "to_bc_feed_tank": {"H2O": 183.5 / 305.1, "TDS": 183.5 / 305.1},
        "to_conc_waste": {"H2O": 6.4 / 305.1, "TDS": 6.4 / 305.1},
    }
    set_separator_op_conditions(m.fs.ro_reject_tank, split_fractions=ro_rej_tank_splits)

    set_pump_op_conditions(m.fs.ro_pump, pressure=400)

    ro_splits = {"to_ro_permeate": {"TDS": 0.01}}
    set_ro_op_conditions(m.fs.ro, split_fractions=ro_splits)

    bc_splits = {}

    split_ab = 350 / (350 * 2 + 450)
    split_c = 1 - 2 * split_ab

    for i, bc_label in enumerate(m.bc_outlets):
        bc_splits[bc_label] = {}
        if i == 2:
            continue
        else:
            bc_splits[bc_label]["H2O"] = split_ab
            bc_splits[bc_label]["TDS"] = split_ab

    set_separator_op_conditions(m.fs.bcs, split_fractions=bc_splits)

    if m.add_basic_bcs:
        splits = {"to_conc_waste": {"H2O": 0.05, "TDS": 0.99}}
        set_bcs_basic_op_conditions(m, splits=splits)


def add_bcs(m, bc_recovery=0.92, ro_recovery=0.7):

    for bc_label in m.BCs:
        bc_fs = FlowsheetBlock()
        m.fs.add_component(bc_label, bc_fs)
        bc_fs = m.fs.find_component(bc_label)
        build_bc(m, bc_fs)
        add_bc_costing(bc_fs)

        # print(f"dof {bc_label} = {degrees_of_freedom(bc_fs)}")
        set_bc_operating_conditions(bc_fs)
        set_bc_scaling(bc_fs)

        feed_props = m.fs.bcs.unit.find_component(f"to_{bc_label.lower()}_state")
        upstream = m.fs.bcs.unit.find_component(f"to_{bc_label.lower()}")
        upstream_to_bc = Arc(source=upstream, destination=bc_fs.feed.inlet)
        init_bc(bc_fs, feed_props=feed_props[0])
        bc_fs.add_component(f"upstream_to_{bc_label}", upstream_to_bc)
        propagate_state(upstream_to_bc)
        TransformationFactory("network.expand_arcs").apply_to(m)
        bc_fs.feed.initialize()
        bc_fs.recovery_mass.unfix()
        bc_fs.recovery_vol.fix(bc_recovery)
        results = solve_bc(bc_fs)

    m.fs.ro_pump.unit.control_volume.properties_out[0].pressure.unfix()
    m.fs.ro.unit.split_fraction[0, "to_ro_permeate", "H2O"].fix(ro_recovery)
    results = solve_bc(m)

    for bc_label in m.BCs:
        bc_fs = m.fs.find_component(bc_label)
        # print(f"dof {bc_label} = {degrees_of_freedom(bc_fs)}")
        bc_fs.compressor.pressure_ratio.fix(1.6)

    # print(f"dof = {degrees_of_freedom(m)}")

    results = solve_bc(m)

    for bc_label in m.BCs:
        bc_fs = m.fs.find_component(bc_label)
        conc_waste_port = m.fs.conc_waste.find_component(f"from_{bc_label.lower()}")
        bc_to_conc_waste = Arc(
            source=bc_fs.disposal.outlet, destination=conc_waste_port.inlet
        )
        m.fs.add_component(f"{bc_label.lower()}_to_conc_waste", bc_to_conc_waste)
        bc_to_conc_waste = m.fs.find_component(f"{bc_label.lower()}_to_conc_waste")
        propagate_state(bc_to_conc_waste)

        demin_port = m.fs.demin.find_component(f"from_{bc_label.lower()}")
        bc_to_demin = Arc(source=bc_fs.product.outlet, destination=demin_port.inlet)
        m.fs.add_component(f"{bc_label.lower()}_to_demin", bc_to_demin)
        bc_to_demin = m.fs.find_component(f"{bc_label.lower()}_to_demin")
        propagate_state(bc_to_demin)

    TransformationFactory("network.expand_arcs").apply_to(m)
    # print(f"dof = {degrees_of_freedom(m)}")

    iscale.calculate_scaling_factors(m)

    for i, bc_label in enumerate(m.BCs):
        if i == 0:
            bc0 = m.fs.find_component(bc_label)
            continue
        bc_fs = m.fs.find_component(bc_label)
        evap_area_constr = Constraint(expr=bc0.evaporator.area == bc_fs.evaporator.area)
        hx_brine_area_constr = Constraint(expr=bc0.hx_brine.area == bc_fs.hx_brine.area)
        hx_distillate_area_constr = Constraint(
            expr=bc0.hx_distillate.area == bc_fs.hx_distillate.area
        )
        m.fs.add_component(f"{bc_label}_evap_area_constr", evap_area_constr)
        m.fs.add_component(f"{bc_label}_hx_brine_area_constr", hx_brine_area_constr)
        m.fs.add_component(
            f"{bc_label}_hx_distillate_area_constr", hx_distillate_area_constr
        )

    # print(f"dof = {degrees_of_freedom(m)}")
    m.fs.costing.cost_process()
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
    )


def add_bcs_basic(m):
    for bc_label in m.BCs:
        m.fs.add_component(bc_label.lower(), FlowsheetBlock(dynamic=False))
        bc_fs = m.fs.find_component(bc_label.lower())
        build_separator(
            bc_fs,
            prop_package=m.fs.properties_feed,
            outlet_list=["to_conc_waste", "to_demin"],
        )
        to_bcs = m.fs.bcs.find_component(f"to_{bc_label.lower()}")
        to_conc_waste = bc_fs.find_component("to_conc_waste")
        to_demin = bc_fs.find_component("to_demin")
        from_demin = m.fs.demin.find_component(f"from_{bc_label.lower()}")
        from_conc_waste = m.fs.conc_waste.find_component(f"from_{bc_label.lower()}")

        a = Arc(source=to_bcs.outlet, destination=bc_fs.feed.inlet)
        m.fs.add_component(f"bcs_to_{bc_label.lower()}", a)

        a = Arc(source=to_demin.outlet, destination=from_demin.inlet)
        m.fs.add_component(f"{bc_label.lower()}_to_demin", a)

        a = Arc(source=to_conc_waste.outlet, destination=from_conc_waste.inlet)
        m.fs.add_component(f"{bc_label.lower()}_to_conc_waste", a)


def set_bcs_basic_op_conditions(m, splits={}):
    if splits == {}:
        raise ValueError("Must provide split fractions for BCs.")
    for bc_label in m.BCs:
        bc_fs = m.fs.find_component(bc_label.lower())
        set_separator_op_conditions(bc_fs, split_fractions=splits)
        print(f"dof {degrees_of_freedom(m)} after setting {bc_label} BCs.")


def init_bcs_basic(m):
    for bc_label in m.BCs:

        bc_fs = m.fs.find_component(bc_label.lower())
        a = m.fs.find_component(f"bcs_to_{bc_label.lower()}")
        propagate_state(a)

        init_separator(bc_fs)
        a = m.fs.find_component(f"{bc_label.lower()}_to_demin")
        propagate_state(a)

        a = m.fs.find_component(f"{bc_label.lower()}_to_conc_waste")
        propagate_state(a)


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
    propagate_state(m.fs.cooling_tower_to_ww_surge_tank)

    init_product(m.fs.cooling_tower_evap)
    init_product(m.fs.service_and_fire)

    init_separator(m.fs.ww_surge_tank)
    propagate_state(m.fs.ww_surge_tank_to_ro_reject_tank)
    propagate_state(m.fs.ww_surge_tank_to_bc_feed_tank)
    propagate_state(m.fs.ww_surge_tank_to_mystery)

    init_product(m.fs.mystery)

    init_separator(m.fs.ro_reject_tank)
    propagate_state(m.fs.ro_reject_tank_to_bc_feed_tank)
    propagate_state(m.fs.ro_reject_tank_to_ro_containment)
    propagate_state(m.fs.ro_reject_tank_to_uf)
    propagate_state(m.fs.ro_reject_tank_to_conc_waste)

    init_separator(m.fs.uf)
    propagate_state(m.fs.uf_to_ro_containment)
    propagate_state(m.fs.uf_to_ro_pump)

    init_pump(m.fs.ro_pump)
    propagate_state(m.fs.ro_pump_to_ro)

    init_ro(m.fs.ro)
    propagate_state(m.fs.ro_to_ro_permeate)
    propagate_state(m.fs.ro_to_ro_containment)

    init_mixer(m.fs.ro_containment)
    propagate_state(m.fs.ro_containment_to_evaporation_ponds)

    init_mixer(m.fs.bc_feed_tank)
    propagate_state(m.fs.bc_feed_tank_to_bcs)

    # if m.add_basic_bcs:
    init_separator(m.fs.bcs)

    init_mixer(m.fs.conc_waste)
    propagate_state(m.fs.conc_waste_to_evaporation_ponds)
    propagate_state(m.fs.ro_containment_to_evaporation_ponds)

    init_mixer(m.fs.evaporation_ponds)
    propagate_state(m.fs.evaporation_ponds_to_pond_evap)

    init_mixer(m.fs.demin)
    propagate_state(m.fs.demin_to_product)

    init_separator(m.fs.permeate)
    propagate_state(m.fs.permeate_to_raw_water_tank_mix)
    propagate_state(m.fs.ro_permeate_to_demin)

    init_product(m.fs.pond_evap)

    if m.add_basic_bcs:
        init_bcs_basic(m)
    # else:
    #     init_bcs(m)

    m.fs.product.initialize()


def print_stream_flows(m, w=30):

    for fb in m.fs.component_objects(Block, descend_into=False):

        if any(x in fb.name for x in ["_expanded", "costing"]):
            continue
        if fb is m.fs.properties_feed:
            continue
        if fb is m.fs.properties_vapor:
            continue
        if "bc_" in fb.name.lower():
            continue

        title = fb.name.replace("fs.", "").replace("_", " ").upper()
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "=" * side + f" {title} " + "=" * side
        print(f"\n{header}\n")
        # print(f"{title.lower()} DOF = {degrees_of_freedom(fb)}\n")

        if (
            isinstance(fb, Feed)
            or isinstance(fb, Product)
            or isinstance(fb, StateJunction)
        ):
            props = fb.find_component("properties")
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
            continue

        u = fb.find_component("unit")

        if u is None:
            continue

        # print(type(u))
        if (
            isinstance(u, Feed)
            or isinstance(u, Product)
            or isinstance(u, StateJunction)
        ):
            props = u.find_component("properties")
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
            continue
        elif isinstance(u, Separator) or isinstance(u, Mixer):
            if u.name == "fs.ro":
                recov = u.split_fraction[0, "to_ro_permeate", "H2O"]() * 100
                print(f'{"Recovery":<{w}s}{f"{recov:<{w},.2f}"}{"%":<{w}s}')
            ms = u.find_component("mixed_state")
            if isinstance(u, Separator):
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
                        pyunits.convert(
                            u.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                            to_units=pyunits.gallons / pyunits.minute,
                        )
                    )
                    for x in u.config.outlet_list
                )
                print(
                    f'{"TOTAL OUTLET FLOW":<{w}s}{f"{tot_flow_out:<{w},.1f}"}{"gpm":<{w}s}'
                )
                for x in u.config.outlet_list:
                    sb = u.find_component(f"{x}_state")
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
                continue
            elif isinstance(u, Mixer):
                tot_flow_in = sum(
                    value(
                        pyunits.convert(
                            u.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                            to_units=pyunits.gallons / pyunits.minute,
                        )
                    )
                    for x in u.config.inlet_list
                )
                print(
                    f'{"TOTAL INLET FLOW":<{w}s}{f"{tot_flow_in:<{w},.1f}"}{"gpm":<{w}s}'
                )
                for x in u.config.inlet_list:
                    sb = u.find_component(f"{x}_state")
                    flow_in = value(
                        pyunits.convert(
                            sb[0].flow_vol_phase["Liq"],
                            to_units=pyunits.gallons / pyunits.minute,
                        )
                    )
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

        elif isinstance(u, Pump):
            cv = u.find_component("control_volume")
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
            print(f'{"Outlet Pressure":<{w}s}{f"{pressure:<{w},.1f}"}{"bar":<{w}s}')
            print(f'{"Outlet Flow":<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}')
            print(f'{"Outlet TDS":<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}')
    for bc_label in m.BCs:
        title = bc_label.replace("fs.", "").replace("_", " ").upper()
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "=" * side + f" {title} " + "=" * side
        print(f"\n{header}\n")
        bc_fs = m.fs.find_component(bc_label)
        print_bc_stream_flows(bc_fs, w=w)


def run_srp_basic():

    m = build_srp(add_basic_bcs=True)
    connect_srp(m)
    set_srp_scaling(m)
    set_srp_operating_conditions(m)
    initialize_srp(m)
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    # print_stream_flows(m)

    return m


def run_srp():
    m = build_srp()
    connect_srp(m)
    set_srp_scaling(m)
    set_srp_operating_conditions(m)
    initialize_srp(m)
    add_bcs(m)
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    print_stream_flows(m)

    return m


if __name__ == "__main__":
    # run_srp()
    run_srp_basic()
