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

from watertap.unit_models.pressure_changer import Pump
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.solvers import get_solver

from srp.utils.utils import touch_flow_and_conc

__all__ = [
    "build_pump",
    "set_pump_op_conditions",
    "init_pump",
    "set_pump_scaling"
]

solver = get_solver()


def build_system(Qin=11343, Cin=1467, feed_temp=27):

    m = ConcreteModel()

    m.Qin = Qin * pyunits.gallons / pyunits.minute
    m.Cin = Cin * pyunits.mg / pyunits.L
    m.feed_temp = feed_temp

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.pump = FlowsheetBlock(dynamic=False)
    build_pump(m.fs.pump, prop_package=m.fs.properties)

    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.feed_to_pump = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.pump.feed.inlet,
    )

    m.fs.pump_to_product = Arc(
        source=m.fs.pump.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_pump(blk, name=None, prop_package=None):

    if name is None:
        name = blk.name.split(".")[-1]

    name = name.replace("_", " ").upper()

    print(f'\n{f"=======> BUILDING {name} UNIT <=======":^60}\n')

    if prop_package is None:
        prop_package = blk.model().fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = Pump(property_package=prop_package)
    touch_flow_and_conc(blk.unit)

    blk.unit.op_pressure_slope = Param(
        initialize=15.7,
        mutable=True,
        units=pyunits.psi / (pyunits.kg / pyunits.m**3),
        doc="Slope of pump outlet pressure vs salinity",
    )
    blk.unit.op_pressure_intercept = Param(
        initialize=96.47,
        mutable=True,
        units=pyunits.psi,
        doc="Intercept of pump outlet pressure vs salinity",
    )
    blk.unit.pressure_constr = Constraint(
        expr=blk.unit.control_volume.properties_out[0].pressure
        == pyunits.convert(
            blk.unit.op_pressure_slope
            * blk.unit.control_volume.properties_in[0].conc_mass_phase_comp[
                "Liq", "TDS"
            ]
            + blk.unit.op_pressure_intercept,
            to_units=pyunits.Pa,
        )
    )
    blk.unit.pressure_constr.deactivate()

    blk.product = StateJunction(property_package=prop_package)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_pump_scaling(blk):
    iscale.set_scaling_factor(blk.unit.work_mechanical, 1e-3)


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

    set_pump_scaling(m.fs.pump)

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


def set_pump_op_conditions(blk, pressure=350, efficiency=0.8):

    pressure = pressure * pyunits.psi

    blk.unit.efficiency_pump.fix(efficiency)
    blk.unit.control_volume.properties_out[0].pressure.fix(pressure)


def init_pump(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()
    propagate_state(blk.unit_to_product)

    blk.product.initialize()


def init_system(m):

    m.fs.feed.initialize()

    propagate_state(m.fs.feed_to_pump)
    init_pump(m.fs.pump)

    propagate_state(m.fs.pump_to_product)
    m.fs.product.initialize()


def main():
    m = build_system()
    set_system_scaling(m)
    set_system_op_conditions(m)
    set_pump_op_conditions(m.fs.pump, pressure=350)
    init_system(m)

    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)


if __name__ == "__main__":
    main()
