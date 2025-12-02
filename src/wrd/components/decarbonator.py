from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Var,
    Objective,
    NonNegativeReals,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
)
from idaes.models.unit_models import StateJunction
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core.util.initialization import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.decarb_system = FlowsheetBlock(dynamic=False)
    build_decarbonator(m.fs.decarb_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_decarbonator(blk, prop_package):

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.unit = StateJunction(property_package=prop_package)

    blk.unit.power_consumption = Var(
        initialize=0,
        domain=NonNegativeReals,
        units=pyunits.kW,
        doc="Power consumption of decarbonator",
    )

    # Add Arcs
    blk.feed_to_decarb = Arc(source=blk.feed.outlet, destination=blk.unit.inlet)
    blk.decarb_to_product = Arc(source=blk.unit.outlet, destination=blk.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_inlet_conditions(blk, Qin=0.581, Cin=0, P_in=10.3):
    """
    Set the operation conditions for the decarbonator
    """
    Qin = (
        Qin * pyunits.m**3 / pyunits.s
    )  # Feed flow rate in m3/s # 0.581 m3/s datapoint from 2/1/2021 0:01. approx. max flow
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    blk.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed.properties[0].pressure.fix(P_in * pyunits.bar)
    blk.feed.properties[0].flow_vol  # Touching

    m = blk.model()
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def set_decarbonator_op_conditions(blk):
    # Instead of hard coding value, this could be in yaml config file
    blk.unit.power_consumption.fix(2.5 * pyunits.kW)
    # Decarb fit params: y = a*x + b
    # blk.unit.power_eq = Constraint(
    #     expr = blk.unit.power_consumption == a * blk.unit.feed.properties[0].flow_vol + b
    # )


def add_decarbonator_scaling(blk):
    set_scaling_factor(blk.unit.power_consumption, 1e-3)


def initialize_decarbonator(blk):
    blk.feed.initialize()
    propagate_state(blk.feed_to_decarb)
    blk.unit.initialize()
    propagate_state(blk.decarb_to_product)
    blk.product.initialize()


def cost_decarbonator(blk):
    blk.costing_package.cost_flow(blk.unit.power_consumption, "electricity")


def report_decarbonator(blk, w=30):
    title = "Decarbonator Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    total_flow = blk.feed.properties[0].flow_vol
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


if __name__ == "__main__":
    m = build_system()  # optional input of stage_num
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_inlet_conditions(m.fs.decarb_system)
    set_decarbonator_op_conditions(m.fs.decarb_system)
    add_decarbonator_scaling(m.fs.decarb_system)
    calculate_scaling_factors(m)
    print(
        f"{degrees_of_freedom(m)} degrees of freedom after setting op and inlet conditions"
    )
    initialize_decarbonator(m.fs.decarb_system)
    report_decarbonator(m.fs.decarb_system)
    # cost_decarbonator(m.fs.decarb_system) # Haven't done any costing yet
    # m.fs.obj = Objective(
    #     expr=m.fs.decarb_system.feed.properties[0].flow_vol_phase["Liq"]
    # )  # There is no D.o.f to optimize with
    # solver = get_solver()
    # results = solver.solve(m)
    # assert_optimal_termination(results)
    
