from pyomo.environ import ConcreteModel, TransformationFactory, value, units as pyunits
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

__all__ = ["build_product", "init_product"]

solver = get_solver()


def build_system(Qin=11343, Cin=1467, split=0.5, feed_temp=27):

    m = ConcreteModel()

    m.Qin = Qin * pyunits.gallons / pyunits.minute
    m.Cin = Cin * pyunits.mg / pyunits.L
    m.feed_temp = feed_temp

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = FlowsheetBlock(dynamic=False)
    build_product(m.fs.product, name="generic_product", prop_package=m.fs.properties)
    m.fs.feed_to_product = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.product.feed.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


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

    iscale.calculate_scaling_factors(m)


def build_product(blk, name=None, prop_package=None):

    if name is None:
        name = blk.name.split(".")[-1]

    name = name.replace("_", " ").upper()

    print(f'\n{f"=======> BUILDING {name} UNIT <=======":^60}\n')

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.product = Product(property_package=prop_package)
    touch_flow_and_conc(blk.product)

    blk.feed_to_product = Arc(
        source=blk.feed.outlet,
        destination=blk.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)

    return blk


def init_product(blk):

    print(f'\n{f"=======> INITIALIZING PRODUCT UNIT <=======":^60}\n')

    propagate_state(blk.feed_to_product)
    blk.product.initialize()


if __name__ == "__main__":
    m = build_system()
    set_system_scaling(m)
    set_system_op_conditions(m)
    print(f"dof = {degrees_of_freedom(m)}")
    init_product(m.fs.product)
    solver = get_solver()
    results = solver.solve(m, tee=True)
