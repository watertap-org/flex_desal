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
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.models.unit_models import Feed, StateJunction, Product

from watertap.unit_models.pressure_changer import Pump
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.solvers import get_solver

from srp.utils.utils import touch_flow_and_conc

__all__ = ["build_pump", "set_pump_op_conditions", "init_pump", "set_pump_scaling"]

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

    # TODO: enable setting of outlet pressure by salinity
    # Not sure if this makes sense they way it is currently implemented
    # blk.unit.op_pressure_slope = Param(
    #     initialize=15.7,
    #     mutable=True,
    #     units=pyunits.psi / (pyunits.kg / pyunits.m**3),
    #     doc="Slope of pump outlet pressure vs salinity",
    # )
    # blk.unit.op_pressure_intercept = Param(
    #     initialize=96.47,
    #     mutable=True,
    #     units=pyunits.psi,
    #     doc="Intercept of pump outlet pressure vs salinity",
    # )
    # blk.unit.pressure_constr = Constraint(
    #     expr=blk.unit.control_volume.properties_out[0].pressure
    #     == pyunits.convert(
    #         blk.unit.op_pressure_slope
    #         * blk.unit.control_volume.properties_in[0].conc_mass_phase_comp[
    #             "Liq", "TDS"
    #         ]
    #         + blk.unit.op_pressure_intercept,
    #         to_units=pyunits.Pa,
    #     ),
    #     doc="Pump outlet pressure constraint based on salinity",
    # )

    # # Deactivated by default to set pressure manually
    # blk.unit.pressure_constr.deactivate()

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


def set_pump_op_conditions(blk, pressure=350, efficiency=0.8, use_pressure_constr=True):
    pressure = pressure * pyunits.psi

    blk.unit.efficiency_pump.fix(efficiency)

    # if use_pressure_constr:
    #     blk.unit.pressure_constr.activate()
    #     blk.unit.control_volume.properties_out[0].pressure.unfix()
    # else:
    blk.unit.control_volume.properties_out[0].pressure.fix(pressure)


def init_pump(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    # if blk.unit.pressure_constr.active:
    #     calculate_variable_from_constraint(
    #         blk.unit.control_volume.properties_out[0].pressure,
    #         blk.unit.pressure_constr,
    #     )
    blk.unit.initialize()
    propagate_state(blk.unit_to_product)

    blk.product.initialize()


def init_system(m):

    m.fs.feed.initialize()

    propagate_state(m.fs.feed_to_pump)
    init_pump(m.fs.pump)

    propagate_state(m.fs.pump_to_product)
    m.fs.product.initialize()


def report_pump(blk, w=35):

    pump_power_watt = value(pyunits.convert(blk.work_mechanical[0], to_units=pyunits.W))

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
    print(f'{"Pump Flow":<{w}s}{value(pump_flow_out):<{w}.1f}{"gpm":<{w}s}')

    print(f'{"Pump Power":<{w}s}{value(pump_power_watt):<{w}.1f}{"W":<{w}s}')
    print(
        f'{"∆P":<{w}s}{value(blk.deltaP[0]):<{w}.1f}{f"{pyunits.get_units(blk.deltaP[0])}":<{w}s}'
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


def main(pressure=350, efficiency=0.8):
    m = build_system()
    set_system_scaling(m)
    set_system_op_conditions(m)
    set_pump_op_conditions(
        m.fs.pump,
        pressure=pressure,
        efficiency=efficiency,
    )
    init_system(m)

    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)

    report_pump(m.fs.pump.unit)

    return m


if __name__ == "__main__":
    m = main()
    m.fs.pump.unit.work_mechanical.display()
    m.fs.pump.unit.control_volume.properties_out[0].pressure.display()
