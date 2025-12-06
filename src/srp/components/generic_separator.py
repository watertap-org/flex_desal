from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.models.unit_models import Separator, StateJunction
from idaes.models.unit_models.separator import SplittingType

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.solvers import get_solver

from srp.utils.utils import touch_flow_and_conc

__all__ = [
    "build_separator",
    "init_separator",
    "set_separator_op_conditions",
    "report_separator",
]


def build_system(Qin=11343, Cin=1467, feed_temp=27, outlet_list=["outlet1", "outlet2"]):

    m = ConcreteModel()

    m.Qin = Qin * pyunits.gallons / pyunits.minute
    m.Cin = Cin * pyunits.mg / pyunits.L
    m.feed_temp = feed_temp

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.separator = FlowsheetBlock(dynamic=False)
    build_separator(
        m.fs.separator, outlet_list=outlet_list, prop_package=m.fs.properties
    )
    m.fs.feed_to_separator = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.separator.feed.inlet,
    )

    for i, outlet in enumerate(outlet_list, 1):
        m.fs.add_component(f"product{i}", Product(property_package=m.fs.properties))
        prod = m.fs.find_component(f"product{i}")
        src = m.fs.separator.find_component(f"{outlet}")
        m.fs.add_component(
            f"{outlet}_to_product{i}", Arc(source=src.outlet, destination=prod.inlet)
        )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_separator(blk, name=None, prop_package=None, outlet_list=None):

    if name is None:
        name = blk.name.split(".")[-1]

    name = name.replace("_", " ").upper()

    print(f'\n{f"=======> BUILDING {name} UNIT <=======":^60}\n')

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    if outlet_list is None:
        raise ValueError(f"Please provide outlet_list for {blk.name}")

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = Separator(
        property_package=prop_package,
        outlet_list=outlet_list,
        split_basis=SplittingType.componentFlow,
    )
    touch_flow_and_conc(blk.unit)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    for outlet in outlet_list:
        # Add StateJunction for each outlet
        blk.add_component(f"{outlet}", StateJunction(property_package=prop_package))
        sj = blk.find_component(f"{outlet}")
        touch_flow_and_conc(sj)

        outlet_port = blk.unit.find_component(f"{outlet}")

        # Connect unit outlet to StateJunction
        blk.add_component(
            f"unit_to_{outlet}",
            Arc(
                source=outlet_port,
                destination=sj.inlet,
            ),
        )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_system_scaling(m):

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s) * 1000),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s)),
        index=("Liq", "TDS"),
    )

    iscale.calculate_scaling_factors(m)


def set_system_op_conditions(m):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): m.Qin,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.Cin,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + m.feed_temp,
        },
        hold_state=True,
    )


def set_separator_op_conditions(blk, split_fractions={}):
    """
    Set split fractions for separator unit.

    split_fractions must take the form
        {
            "outlet1": {"H2O": 0.3, "TDS": 0.3},
            "outlet2": {"H2O": 0.1, "TDS": 0.1},
        }
    There should be one split fraction that is left unspecified to avoid over-specification.
    """

    if split_fractions == {}:
        raise ValueError(f"Please provide split_fractions dict for {blk.name} unit")

    for port, comp_splits in split_fractions.items():
        for comp, frac in comp_splits.items():
            if frac is not None:
                blk.unit.split_fraction[0, port, comp].fix(frac)


def init_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_separator)
    init_separator(m.fs.separator)

    for i, outlet in enumerate(m.fs.separator.unit.config.outlet_list, 1):
        a = m.fs.find_component(f"{outlet}_to_product{i}")
        propagate_state(a)
        pb = m.fs.find_component(f"product{i}")
        pb.initialize()


def init_separator(blk, name=None):

    if name is None:
        name = blk.name.split(".")[-1]

    name = name.replace("_", " ").upper()

    print(f'\n{f"=======> INITIALIZING {name} UNIT <=======":^60}\n')

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)
    blk.unit.initialize()

    for outlet in blk.unit.config.outlet_list:
        a = blk.find_component(f"unit_to_{outlet}")
        propagate_state(a)
        pb = blk.find_component(f"{outlet}")
        pb.initialize()


def report_separator(blk, w=25):

    fs = blk.flowsheet()
    title = fs.name.replace("fs.", "").replace("_", " ").upper()

    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    ms = blk.find_component("mixed_state")
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
    print(f'{"INLET Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"INLET TDS":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
    tot_flow_out = sum(
        value(
            value(
                pyunits.convert(
                    blk.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                    to_units=pyunits.gallons / pyunits.minute,
                )
            )
        )
        for x in blk.config.outlet_list
    )
    print(f'{"TOTAL OUTLET FLOW":<{w}s}{f"{tot_flow_out:<{w},.1f}"}{"gpm":<{w}s}')
    for x in blk.config.outlet_list:
        sb = blk.find_component(f"{x}_state")
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


def main():

    m = build_system(outlet_list=["cats", "dogs"])
    set_system_scaling(m)
    set_system_op_conditions(m)

    split_fractions = {"cats": {"H2O": 0.9, "TDS": 0.9}}

    set_separator_op_conditions(
        m.fs.separator,
        split_fractions=split_fractions,
    )
    init_system(m)
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    print(f"dof = {degrees_of_freedom(m)}")
    report_separator(m.fs.separator.unit)


if __name__ == "__main__":
    main()
