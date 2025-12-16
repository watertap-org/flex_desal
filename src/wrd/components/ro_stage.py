from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    units as pyunits,
    TransformationFactory,
)
from pyomo.network import Arc

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed, Product
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import StateJunction
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.costing import WaterTAPCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from wrd.utilities import load_config, get_config_file
from wrd.components.pump import *
from wrd.components.ro import *
from srp.utils import touch_flow_and_conc

__all__ = [
    "build_ro_stage",
    "set_ro_stage_op_conditions",
    "set_ro_stage_scaling",
    "initialize_ro_stage",
    "report_ro_stage",
    "add_ro_stage_costing",
]

solver = get_solver()

"""
Component file for single RO stage consisting of a pump and RO unit.

feed > pump > RO > product
                 > disposal
"""


def build_system(stage_num=1, file="wrd_ro_inputs_8_19_21.yaml"):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.ro_stage = FlowsheetBlock(dynamic=False)
    build_ro_stage(
        m.fs.ro_stage, stage_num=stage_num, file=file, prop_package=m.fs.properties
    )

    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.brine = Product(property_package=m.fs.properties)

    # Arcs to connect the unit models
    m.fs.feed_to_ro_stage = Arc(
        source=m.fs.feed.outlet, destination=m.fs.ro_stage.feed.inlet
    )
    m.fs.ro_stage_to_product = Arc(
        source=m.fs.ro_stage.product.outlet, destination=m.fs.product.inlet
    )
    m.fs.ro_stage_to_brine = Arc(
        source=m.fs.ro_stage.disposal.outlet, destination=m.fs.brine.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def set_inlet_conditions(m, Qin=2637, Cin=0.5, Tin=298, Pin=101325):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): Pin,
            ("temperature", None): Tin,
        },
        hold_state=True,
    )


def set_ro_stage_op_conditions(blk):

    set_pump_op_conditions(blk.pump)
    set_ro_op_conditions(blk.ro)


def set_ro_stage_scaling(blk):

    add_pump_scaling(blk.pump)
    set_ro_scaling(blk.ro)


def build_ro_stage(
    blk,
    stage_num=1,
    file="wrd_ro_inputs_8_19_21.yaml",
    prop_package=None,
):
    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.config_data = load_config(get_config_file(file))

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.pump = FlowsheetBlock(dynamic=False)
    build_pump(blk.pump, stage_num=stage_num, file=file, prop_package=prop_package)

    blk.ro = FlowsheetBlock(dynamic=False)
    build_ro(blk.ro, stage_num=stage_num, file=file, prop_package=prop_package)

    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    # Arcs to connect the unit models
    blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump.feed.inlet)
    blk.pump_to_ro = Arc(source=blk.pump.product.outlet, destination=blk.ro.feed.inlet)
    blk.ro_to_product = Arc(source=blk.ro.product.outlet, destination=blk.product.inlet)
    blk.ro_to_disposal = Arc(
        source=blk.ro.disposal.outlet, destination=blk.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def initialize_system(m):
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_ro_stage)

    initialize_ro_stage(m.fs.ro_stage)

    propagate_state(m.fs.ro_stage_to_product)
    m.fs.product.initialize()
    propagate_state(m.fs.ro_stage_to_brine)
    m.fs.brine.initialize()


def initialize_ro_stage(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_pump)

    initialize_pump(blk.pump)

    propagate_state(blk.pump_to_ro)
    initialize_ro(blk.ro)

    propagate_state(blk.ro_to_product)
    blk.product.initialize()
    propagate_state(blk.ro_to_disposal)
    blk.disposal.initialize()


def report_ro_stage(blk, w=30):
    # title = "RO Stage Report"
    # side = int(((3 * w) - len(title)) / 2) - 1
    # header = "=" * side + f" {title} " + "=" * side
    # print(f"\n{header}\n")
    report_pump(blk.pump, w=w)
    report_ro(blk.ro, w=w)


def add_ro_stage_costing(blk, costing_package=None):

    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    add_pump_costing(blk.pump, costing_package=costing_package)
    add_ro_costing(blk.ro, costing_package=costing_package)


def main(
    Qin=2637,
    Cin=0.528,
    Tin=302,
    Pin=101325,
    stage_num=1,
    file="wrd_ro_inputs_8_19_21.yaml",
):

    m = build_system(stage_num=stage_num, file=file)
    set_ro_stage_scaling(m.fs.ro_stage)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin, Tin=Tin, Pin=Pin)
    set_ro_stage_op_conditions(m.fs.ro_stage)

    add_ro_stage_costing(m.fs.ro_stage)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
    )

    initialize_system(m)

    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_ro_stage(m.fs.ro_stage)

    return m


def run_august_stages():
    # August 19, 2021 Data
    # Stage 1
    m = main()
    # Stage 2
    m = main(
        Qin=1029,
        Cin=1.2479,
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_ro_inputs_8_19_21.yaml",
    )
    # Stage 3
    m = main(
        Qin=384,
        Cin=4.847 / 2,
        Tin=302,
        Pin=(112.6 - 41.9) * pyunits.psi,
        stage_num=3,
        file="wrd_ro_inputs_8_19_21.yaml",
    )


def run_march_stages():

    # March 13, 2021 Data
    # Stage 1
    m = main(
        Qin=2452,
        Cin=0.503,
        Tin=295,
        Pin=101325,
        stage_num=1,
        file="wrd_ro_inputs_3_13_21.yaml",
    )
    # Stage 2
    m = main(
        Qin=1047,
        Cin=2.239 / 2,
        Tin=295,
        Pin=143.5 * pyunits.psi,
        stage_num=2,
        file="wrd_ro_inputs_3_13_21.yaml",
    )
    # Stage 3
    m = main(
        Qin=506.5,
        Cin=2.7,
        Tin=295,
        Pin=(106.3 - 59.3) * pyunits.psi,
        stage_num=3,
        file="wrd_ro_inputs_3_13_21.yaml",
    )


if __name__ == "__main__":
    run_august_stages()
    run_march_stages()
