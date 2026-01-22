from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import StateJunction, Feed
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.costing import WaterTAPCosting
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.unit_models.deep_well_injection import (
    DeepWellInjection as BrineDisposal,
)  # just for fun

from wrd.utilities import load_config, get_config_file, get_config_value
from srp.utils import touch_flow_and_conc

solver = get_solver()

__all__ = [
    "build_brine_disposal",
    "initialize_brine_disposal",
    "add_brine_disposal_costing",
    "report_brine_disposal",
]


def build_system(file="wrd_inputs_8_19_21.yaml"):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2021

    config = get_config_file(file)
    m.fs.config_data = load_config(config)

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.brine_disposal = FlowsheetBlock(dynamic=False)
    build_brine_disposal(m.fs.brine_disposal, file=file, prop_package=m.fs.properties)

    # Arcs to connect the unit models
    m.fs.feed_to_brine_disposal = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.brine_disposal.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def build_brine_disposal(blk, file="wrd_inputs_8_19_21.yaml", prop_package=None):

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.ro_properties

    blk.config_data = load_config(get_config_file(file))

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = BrineDisposal(
        property_package=prop_package,
    )

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)

    return blk


def initialize_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_brine_disposal)

    initialize_brine_disposal(m.fs.brine_disposal)


def initialize_brine_disposal(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()


def set_inlet_conditions(m, Qin=2637, Cin=0.5, Tin=302, Pin=101325):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): Pin,
            ("temperature", None): Tin,
        },
        hold_state=True,
    )


def add_brine_disposal_costing(blk, costing_package=None):
    m = blk.model()
    if costing_package is None:
        costing_package = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=costing_package,
        costing_method_arguments={"cost_method": "as_opex"},
    )

    brine_disposal_cost = get_config_value(
        m.fs.config_data, "brine_disposal_cost", "brine_disposal_cost"
    )
    costing_package.deep_well_injection.dwi_lcow.fix(brine_disposal_cost)


def report_brine_disposal(blk, w=25):
    title = "Brine Disposal Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    flow_in = blk.feed.properties[0].flow_vol_phase["Liq"]
    conc_in = blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
    pressure_in = blk.feed.properties[0].pressure
    temperature_in = blk.feed.properties[0].temperature

    print(
        f'{f"Flow Rate":<{w}s}{value(pyunits.convert(flow_in, to_units=pyunits.gal / pyunits.min)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Concentration":<{w}s}{value(pyunits.convert(conc_in, to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
    )
    print(
        f'{f"Pressure":<{w}s}{value(pyunits.convert(pressure_in, to_units=pyunits.psi)):<{w}.3f}{"psi"}'
    )
    print(
        f'{f"Temperature":<{w}s}{value(pyunits.convert(temperature_in, to_units=pyunits.degK)):<{w}.3f}{"°K"}'
    )

    if blk.unit.find_component("costing") is not None:
        cb = blk.model().fs.costing
        dwi_cost = cb.deep_well_injection.dwi_lcow
        print(
            f'{f"Brine Disposal Unit Cost":<{w}s}{value(pyunits.convert(dwi_cost, to_units=cb.base_currency / pyunits.m**3)):<{w}.3f}{f"{pyunits.get_units(dwi_cost)}":<{w}s}'
        )
        print(
            f'{f"Brine Disposal Opex":<{w}s}{value(blk.unit.costing.variable_operating_cost):<{w}.2f}{f"{pyunits.get_units(blk.unit.costing.variable_operating_cost)}":<{w}s}'
        )


def main(
    Qin=717.7,
    Cin=7.558,
    Tin=302,
    Pin=101325,
    file="wrd_inputs_8_19_21.yaml",
):

    m = build_system(file=file)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin, Tin=Tin, Pin=Pin)
    add_brine_disposal_costing(m.fs.brine_disposal)
    initialize_system(m)
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)

    report_brine_disposal(m.fs.brine_disposal, w=40)

    return m


if __name__ == "__main__":
    m = main()
