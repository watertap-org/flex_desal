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
from idaes.models.unit_models import Feed, Mixer, Separator, StateJunction, Product
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.separator import SplittingType

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.solvers import get_solver

from srp.utils.utils import touch_flow_and_conc

solver = get_solver()

__all__ = [
    "build_mixer",
    "init_mixer",
]


def build_system(
    Qin=11343,
    Cin=1467,
    feed_temp=27,
    inlet_list=["inlet1", "inlet2"],
):

    m = ConcreteModel()
    m.inlet_list = inlet_list

    m.Qin = Qin * pyunits.gallons / pyunits.minute
    m.Cin = Cin * pyunits.mg / pyunits.L
    m.feed_temp = feed_temp

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.mixer = FlowsheetBlock(dynamic=False)

    build_mixer(m.fs.mixer, inlet_list=inlet_list)

    for inlet in inlet_list:
        m.fs.add_component(f"{inlet}", Feed(property_package=m.fs.properties))
        feed_unit = m.fs.find_component(f"{inlet}")
        touch_flow_and_conc(feed_unit)
        inlet_port = m.fs.mixer.unit.find_component(f"{inlet}")
        m.fs.add_component(
            f"{inlet}_to_mixer",
            Arc(
                source=feed_unit.outlet,
                destination=inlet_port,
            ),
        )

    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.unit_to_product = Arc(
        source=m.fs.mixer.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_mixer(blk, name=None, prop_package=None, inlet_list=["inlet1", "inlet2"]):

    if name is None:
        name = blk.name.split(".")[-1]

    name = name.replace("_", " ").upper()

    print(f'\n{f"=======> BUILDING {name} UNIT <=======":^60}\n')

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.unit = Mixer(
        property_package=prop_package,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=inlet_list,
    )
    touch_flow_and_conc(blk.unit)

    blk.unit.outlet.pressure[0].fix(101325)

    for inlet in inlet_list:
        blk.add_component(
            f"{inlet}",
            StateJunction(property_package=prop_package),
        )
        sj = blk.find_component(f"{inlet}")
        touch_flow_and_conc(sj)
        blk.add_component(
            f"{inlet}_to_unit",
            Arc(
                source=sj.outlet,
                destination=blk.unit.find_component(f"{inlet}"),
            ),
        )

    blk.product = StateJunction(property_package=prop_package)
    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
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

    for inlet in m.inlet_list:
        feed_unit = m.fs.find_component(f"{inlet}")
        feed_unit.properties.calculate_state(
            var_args={
                ("flow_vol_phase", ("Liq")): m.Qin / len(m.inlet_list),
                ("conc_mass_phase_comp", ("Liq", "TDS")): m.Cin,
                ("pressure", None): 101325,
                ("temperature", None): 273.15 + m.feed_temp,
            },
            hold_state=True,
        )


def init_system(m):

    for inlet in m.inlet_list:
        feed_unit = m.fs.find_component(f"{inlet}")
        feed_unit.initialize()
        a = m.fs.find_component(f"{inlet}_to_mixer")
        propagate_state(a)

    init_mixer(m.fs.mixer)

    propagate_state(m.fs.unit_to_product)

    m.fs.product.initialize()


def init_mixer(blk, name=None):
    if name is None:
        name = blk.name.split(".")[-1]

    name = name.replace("_", " ").upper()

    print(f'\n{f"=======> INITIALIZING {name} UNIT <=======":^60}\n')

    for inlet in blk.unit.config.inlet_list:
        sj = blk.find_component(f"{inlet}")
        sj.initialize()
        a = blk.find_component(f"{inlet}_to_unit")
        propagate_state(a)

    blk.unit.initialize()

    propagate_state(blk.unit_to_product)
    blk.product.initialize()


def main():

    m = build_system()
    set_system_scaling(m)
    set_system_op_conditions(m)
    init_system(m)
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)


if __name__ == "__main__":
    main()
