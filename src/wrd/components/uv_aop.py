from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Param,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
)
from idaes.models.unit_models import StateJunction, Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.solvers import get_solver
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.costing import WaterTAPCosting

from srp.utils import touch_flow_and_conc


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.feed = Feed(property_package=m.fs.properties)

    m.fs.uv_aop_system = FlowsheetBlock(dynamic=False)
    build_uv_aop(m.fs.uv_aop_system, prop_package=m.fs.properties, **kwargs)

    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.feed_to_unit = Arc(
        source=m.fs.feed.outlet, destination=m.fs.uv_aop_system.feed.inlet
    )
    m.fs.unit_to_product = Arc(
        source=m.fs.uv_aop_system.unit.outlet, destination=m.fs.product.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    return m


def build_uv_aop(blk, prop_package=None):

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)
    blk.product = StateJunction(property_package=prop_package)
    blk.unit = StateJunction(property_package=prop_package)

    # Look to RO_rework or UV_Surrogate for loading and incorperating the surrogate model here

    # Using all Feb 2021 data
    # UV1: a = 14.78, b = 8.71. If b = 0, best a = 15.86
    # UV2: a = 14.78, b = 7.39. If b = 0, best a = 15.9

    blk.unit.power_eq_slope = Param(
        initialize=15.9,
        mutable=True,
        units=pyunits.kW / (pyunits.Mgallons / pyunits.day),
        doc="Slope of power consumption equation",
    )
    blk.unit.power_eq_intercept = Param(
        initialize=0,
        mutable=True,
        units=pyunits.kW,
        doc="Intercept of power consumption equation",
    )

    blk.unit.power_consumption = Var(
        initialize=0,
        domain=NonNegativeReals,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Power consumption of UV/AOP unit",
    )

    blk.unit.eq_power = Constraint(
        expr=blk.unit.power_consumption
        == pyunits.convert(
            blk.unit.power_eq_slope * blk.feed.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.kW,
        )
        + blk.unit.power_eq_intercept,
        doc="Linear fit to data",
    )

    blk.feed_to_unit = Arc(source=blk.feed.outlet, destination=blk.unit.inlet)
    blk.unit_to_product = Arc(source=blk.unit.outlet, destination=blk.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_uv_aop_op_conditions(blk):
    pass


def set_inlet_conditions(m, Qin=0.27, Cin=0, P_in=1):
    """
    Set the inlet conditions for the UV/AOP system.
    """

    # Feed flow rate in m3/s
    # Flow into one UV unit, near max flow
    Qin = Qin * pyunits.m**3 / pyunits.s
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    m.fs.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    m.fs.feed.properties[0].pressure.fix(P_in * pyunits.bar)


def initialize_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_unit)

    initialize_uv_aop(m.fs.uv_aop_system)

    propagate_state(m.fs.unit_to_product)
    m.fs.product.initialize()


def initialize_uv_aop(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()
    propagate_state(blk.unit_to_product)

    blk.product.initialize()


def add_uv_aop_scaling(blk):
    set_scaling_factor(blk.unit.power_consumption, 1e-3)


def cost_uv_aop(blk, costing_package=None):
    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    # Using this method to cost electricity consumption because it is a ZO model
    costing_package.cost_flow(
        pyunits.convert(blk.unit.power_consumption, to_units=pyunits.kW), "electricity"
    )


def report_uv(blk, w=30):
    title = "UV Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    total_flow = blk.feed.properties[0].flow_vol_phase["Liq"]
    power = blk.unit.power_consumption
    print(
        f'{f"Total Flow Rate (MGD)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.Mgallons /pyunits.day)):<{w}.3f}{"MGD"}'
    )
    print(f'{f"Total Flow Rate (m3/s)":<{w}s}{value(total_flow):<{w}.3e}{"m3/s"}')
    print(
        f'{f"Total Flow Rate (gpm)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Power Consumption (kW)":<{w}s}{value(pyunits.convert(power, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    m = blk.model()
    SEC = m.fs.costing.SEC
    print(
        f'{f"Specific Energy (SEC)":<{w}s}{value(pyunits.convert(SEC, to_units=pyunits.kWh / pyunits.m**3)):<{w}.3f}{"kWh/m3"}'
    )


def main():

    m = build_system()
    set_inlet_conditions(m)
    set_uv_aop_op_conditions(m.fs.uv_aop_system)
    add_uv_aop_scaling(m.fs.uv_aop_system)
    calculate_scaling_factors(m)
    initialize_system(m)
    cost_uv_aop(m.fs.uv_aop_system)
    m.fs.costing.cost_process()
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"],
        name="SEC",
    )
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_uv(m.fs.uv_aop_system, w=40)
    return m


if __name__ == "__main__":
    m = main()
