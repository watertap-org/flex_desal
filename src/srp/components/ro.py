from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Param,
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
from idaes.models.unit_models import Feed, Separator, StateJunction, Product
from idaes.models.unit_models.separator import SplittingType

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.solvers import get_solver

from srp.utils.utils import touch_flow_and_conc

__all__ = [
    "build_ro",
    "init_ro",
    "set_ro_op_conditions",
]

solver = get_solver()


def build_system(
    Qin=11343,
    Cin=1467,
    split=0.5,
    feed_temp=27,
    outlet_list=["to_ro_containment", "to_ro_permeate"],
):

    m = ConcreteModel()

    m.Qin = Qin * pyunits.gallons / pyunits.minute
    m.Cin = Cin * pyunits.mg / pyunits.L
    m.feed_temp = feed_temp

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.ro = FlowsheetBlock(dynamic=False)

    build_ro(m.fs.ro)
    m.fs.feed_to_ro = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.ro.feed.inlet,
    )
    for i, outlet in enumerate(outlet_list, 1):
        m.fs.add_component(f"product{i}", Product(property_package=m.fs.properties))
        prod = m.fs.find_component(f"product{i}")
        src = m.fs.ro.find_component(f"{outlet}")
        m.fs.add_component(
            f"{outlet}_to_product{i}", Arc(source=src.outlet, destination=prod.inlet)
        )
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def build_ro(
    blk,
    prop_package=None,
    pump=None,
    name=None,
    outlet_list=["to_ro_containment", "to_ro_permeate"],
):
    if name is None:
        name = blk.name.split(".")[-1]

    name = name.replace("_", " ").upper()

    print(f'\n{f"=======> BUILDING {name} UNIT <=======":^60}\n')

    # if pump is None:
    #     raise ValueError("A pump unit must be provided to build the RO system.")

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = Separator(
        property_package=prop_package,
        outlet_list=outlet_list,
        split_basis=SplittingType.componentFlow,
    )
    touch_flow_and_conc(blk.unit)

    if pump is not None:

        blk.feed_osm_pressure = Param(
            initialize=9.76,
            units=pyunits.bar,
            mutable=True,
            doc="Feed osmotic pressure",  # estimated from vant Hoff equation for NaCl @ 25C
        )

        blk.unit.recov_base = Param(initialize=0.1370946, mutable=True)
        blk.unit.recov_exp = Param(initialize=0.58863718, mutable=True)

        blk.unit.recovery_constr = Constraint(
            expr=blk.unit.split_fraction[0, "to_ro_permeate", "H2O"]
            == blk.unit.recov_base
            * (
                pyunits.convert(
                    pump.control_volume.properties_out[0].pressure,
                    to_units=pyunits.bar,
                )
                - blk.feed_osm_pressure
            )
            ** blk.unit.recov_exp
        )

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


def set_ro_op_conditions(blk, split_fractions={}):

    blk.unit.split_fraction[0, "to_ro_permeate", "H2O"].setlb(0.5)
    blk.unit.split_fraction[0, "to_ro_permeate", "H2O"].setub(0.8)

    for n, (port, comp_splits) in enumerate(split_fractions.items(), 1):
        for comp, frac in comp_splits.items():
            if frac is not None:
                blk.unit.split_fraction[0, port, comp].fix(frac)


def init_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_ro)

    init_ro(m.fs.ro)

    for i, outlet in enumerate(m.fs.ro.unit.config.outlet_list, 1):
        a = m.fs.find_component(f"{outlet}_to_product{i}")
        propagate_state(a)
        pb = m.fs.find_component(f"product{i}")
        pb.initialize()


def init_ro(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()

    for outlet in blk.unit.config.outlet_list:
        a = blk.find_component(f"unit_to_{outlet}")
        propagate_state(a)
        pb = blk.find_component(f"{outlet}")
        pb.initialize()


def main():
    m = build_system()
    set_system_scaling(m)
    set_system_op_conditions(m)
    split_fractions = {"to_ro_permeate": {"H2O": 0.5, "TDS": 0.5}}
    set_ro_op_conditions(m.fs.ro, split_fractions=split_fractions)
    init_system(m)
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)


if __name__ == "__main__":
    main()
