from pyomo.environ import (
    ConcreteModel,
    Var,
    Param,
    Constraint,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import StateJunction
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
    get_scaling_factor,
)

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.unit_models.pressure_changer import Pump
from watertap.core.solvers import get_solver

from wrd.utilities import load_config, get_config_value, get_config_file
from srp.utils import touch_flow_and_conc

__all__ = [
    "build_pump",
    "initialize_pump",
    "set_pump_op_conditions",
    "report_pump",
    "add_pump_scaling",
]


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.pump_system = FlowsheetBlock(dynamic=False)
    build_pump(m.fs.pump_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_pump(blk, stage_num=1, file="wrd_ro_inputs_8_19_21.yaml", prop_package=None):
    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.ro_properties

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)
    blk.product = StateJunction(property_package=prop_package)
    config_file_name = get_config_file(file)
    blk.config_data = load_config(config_file_name)
    blk.stage_num = stage_num

    blk.unit = Pump(property_package=prop_package)

    # Create variable for the efficiency from the pump curves
    blk.unit.efficiency_fluid = Var(
        initialize=0.7,
        units=pyunits.dimensionless,
        bounds=(0, 1),
        doc="Efficiency from pump curves",
    )

    # Load Values for surrogate model
    if stage_num == 1:
        a_0 = 0.389
        a_1 = -0.535
        a_2 = 41.373
        a_3 = -138.82
    elif stage_num == 2:
        a_0 = 0.067
        a_1 = 21.112
        a_2 = -133.157
        a_3 = -234.386
    else:
        a_0 = 0.067
        a_1 = 21.112
        a_2 = -133.157
        a_3 = -234.386

    # Create Variables for simple "surrogate"
    blk.unit.efficiency_eq_constant = Param(
        initialize=a_0,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Constant term of Efficiency equation",
    )

    blk.unit.efficiency_eq_linear = Param(
        initialize=a_1,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -1,
        doc="Linear term of Efficiency equation",
    )

    blk.unit.efficiency_eq_squared = Param(
        initialize=a_2,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -2,
        doc="Squared term of Efficiency equation",
    )

    blk.unit.efficiency_eq_cubed = Param(
        initialize=a_3,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -3,
        doc="Cubed term of Efficiency equation",
    )

    flow = blk.feed.properties[0].flow_vol_phase["Liq"]

    blk.unit.efficiency_surr_eq = Constraint(
        expr=blk.unit.efficiency_fluid
        == blk.unit.efficiency_eq_cubed * flow**3
        + blk.unit.efficiency_eq_squared * flow**2
        + blk.unit.efficiency_eq_linear * flow
        + blk.unit.efficiency_eq_constant,
        doc="Efficiency surrogate equation",
    )
    blk.unit.efficiency_pump.bounds = (0, 1)

    blk.unit.efficiency_motor = Param(
        initialize=0.938,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Efficiency of motor and VFD",
    )

    blk.unit.efficiency_vfd = Param(
        initialize=0.95,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Efficiency of VFD",
    )

    blk.unit.efficiency_electrical = Constraint(
        expr=blk.unit.efficiency_pump[0]
        == blk.unit.efficiency_motor
        * blk.unit.efficiency_vfd
        * blk.unit.efficiency_fluid
    )

    # Add Arcs
    blk.feed_to_unit = Arc(source=blk.feed.outlet, destination=blk.unit.inlet)
    blk.unit_to_product = Arc(source=blk.unit.outlet, destination=blk.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_pump_op_conditions(blk):
    Pout = get_config_value(
        blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{blk.stage_num}"
    )

    print(f"Setting pump {blk.stage_num} operating conditions, Pout={Pout}")
    blk.unit.control_volume.properties_out[0].pressure.fix(Pout)


def set_inlet_conditions(blk):
    """
    Set the inlet conditions for the Pump
    """
    Qin = get_config_value(
        blk.config_data,
        "pump_flowrate",
        "pumps",
        f"pump_{blk.stage_num}",
    )

    Cin = get_config_value(
        blk.config_data,
        "feed_conductivity",
        "pumps",
        f"pump_{blk.stage_num}",
    ) * get_config_value(blk.config_data, "feed_conductivity_conversion", "feed_stream")

    Pin = get_config_value(
        blk.config_data,
        "pump_suction_pressure",
        "pumps",
        f"pump_{blk.stage_num}",
    )
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    blk.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed.properties[0].pressure.fix(Pin)
    # blk.feed.properties[0].flow_vol  # Touching
    blk.unit.control_volume.properties_in[0].flow_vol  # Touching

    # Scaling properties
    m = blk.model()
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def add_pump_scaling(blk):
    set_scaling_factor(blk.unit.work_mechanical[0], 1e-3)


def initialize_pump(blk):
    # Touch Properties that are needed for later
    blk.feed.properties[0].flow_vol_phase["Liq"]
    blk.product.properties[0].flow_vol_phase["Liq"]

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()
    propagate_state(blk.unit_to_product)
    blk.product.initialize()


def report_pump(blk, w=30):
    title = "Pump Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    flow_in = blk.feed.properties[0].flow_vol_phase["Liq"]
    deltaP = blk.unit.deltaP[0]
    work = blk.unit.work_mechanical[0]
    print(
        f'{f"Inlet Flow":<{w}s}{value(pyunits.convert(flow_in, to_units=pyunits.gallons /pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(f'{f"∆P (Pa)":<{w}s}{value(deltaP):<{w}.3e}{"Pa"}')
    print(
        f'{f"∆P (psi)":<{w}s}{value(pyunits.convert(deltaP, to_units=pyunits.psi)):<{w}.3e}{"psi"}'
    )
    print(
        f'{f"Work Mech. (kW)":<{w}s}{value(pyunits.convert(work, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    print(f'{f"Efficiency (-)":<{w}s}{value(blk.unit.efficiency_pump[0]):<{w}.3f}{"-"}')


def main(stage_num=1, date="8_19_21"):
    m = build_system(stage_num=stage_num, date=date)  # optional input of stage_num
    set_inlet_conditions(m.fs.pump_system)
    set_pump_op_conditions(m.fs.pump_system)
    add_pump_scaling(m.fs.pump_system)
    calculate_scaling_factors(m)
    initialize_pump(m.fs.pump_system)
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    return m


if __name__ == "__main__":
    stage_num = 2
    m = build_system(stage_num=stage_num)  # optional input of stage_num
    assert_units_consistent(m)
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_inlet_conditions(m.fs.pump_system)
    set_pump_op_conditions(m.fs.pump_system)
    print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    add_pump_scaling(m.fs.pump_system)
    calculate_scaling_factors(m)
    initialize_pump(m.fs.pump_system)
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    # print(f"{iscale.jacobian_cond(m.fs.pump_system):.2e}")
    report_pump(m.fs.pump_system)
