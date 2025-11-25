from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
from idaes.models.unit_models import StateJunction
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *

def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.uv_aop_system = FlowsheetBlock(dynamic=False)
    build_uv_aop(m.fs.uv_aop_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_uv_aop(blk, prop_package):
    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.unit = StateJunction(property_package=prop_package)

    blk.unit.power_consumption = Var(
        initialize=0,
        domain=NonNegativeReals,
        units=pyunits.kW,
        doc="Power consumption of UV_aop",
    )
   
    blk.feed_to_UV = Arc(source = blk.feed.inlet,destination= blk.unit.inlet)
    blk.UV_to_product = Arc(source=blk.unit.outlet, destination = blk.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_uv_aop_op_conditions(blk):
    # Look to RO_rework or UV_Surrogate for loading and incorperating the surrogate model here

    # Using all Feb 2021 data
    # UV1: a = 14.78, b = 8.71. If b = 0, best a = 15.86
    # UV2: a = 14.78, b = 7.39. If b = 0, best a = 15.98
    a = 15.9 * pyunits.kW / (pyunits.Mgallons / pyunits.day)
    b = 0 * pyunits.kW

    blk.unit.eq_power = Constraint(
        expr = blk.unit.power_consumption == a*blk.feed.properties[0].flow_vol_phase['Liq'] + b,
        doc = "Linear fit to data"
    )
    

def set_inlet_conditions(blk, Qin= 0.27, Cin= 0, P_in = 10.6):
    """
    Set the operation conditions for the UV. 
    """
    Qin = (Qin) * pyunits.m**3 / pyunits.s  # Feed flow rate in m3/s # Flow into one UV unit, near max flow
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        feed_mass_flow_water
    )
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        feed_mass_flow_salt
    )
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


def initialize_uv_aop(blk):
    blk.feed.initialize()
    propagate_state(blk.feed_to_UV)
    blk.unit.initialize()
    propagate_state(blk.UV_to_product)
    blk.product.initialize()
  
def add_uv_aop_scaling(blk):
    set_scaling_factor(blk.unit.power_consumption,1e-3)

def cost_uv_aop(blk):
    blk.costing_package.cost_flow(blk.unit.power_consumption, "electricity")

def report_uv(blk,w=30):
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

if __name__ == "__main__":
    m = build_system()
    set_inlet_conditions(m.fs.uv_aop_system)
    set_uv_aop_op_conditions(m.fs.uv_aop_system)
    add_uv_aop_scaling(m.fs.uv_aop_system)
    calculate_scaling_factors(m)
    initialize_uv_aop(m.fs.uv_aop_system)
    m.fs.obj = Objective(
        expr=m.fs.uv_aop_system.feed.properties[0].flow_vol_phase["Liq"]
    )
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_uv(m.fs.uv_aop_system, w=40)