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


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.pump_system = FlowsheetBlock(dynamic=False)
    build_wrd_pump(m.fs.pump_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_wrd_pump(blk, stage_num=1, date="8_19_21", prop_package=None):
    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.ro_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    config_file_name = get_config_file("wrd_ro_inputs_" + date + ".yaml")
    blk.config_data = load_config(config_file_name)
    blk.stage_num = stage_num

    blk.pump = Pump(property_package=prop_package)

    # Create variable for the efficiency from the pump curves
    blk.pump.efficiency_fluid = Var(
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
        # Don't have these pump curves yet, so estimating to pass tests
        a_0 = 0.6
        a_1 = 0
        a_2 = 0
        a_3 = 0

    # Create Variables for simple "surrogate"
    blk.pump.efficiency_eq_constant = Param(
        initialize=a_0,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Constant term of Efficiency equation",
    )

    blk.pump.efficiency_eq_linear = Param(
        initialize=a_1,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -1,
        doc="Linear term of Efficiency equation",
    )

    blk.pump.efficiency_eq_squared = Param(
        initialize=a_2,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -2,
        doc="Squared term of Efficiency equation",
    )

    blk.pump.efficiency_eq_cubed = Param(
        initialize=a_3,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -3,
        doc="Cubed term of Efficiency equation",
    )

    flow = blk.feed.properties[0].flow_vol_phase["Liq"]

    blk.pump.efficiency_surr_eq = Constraint(
        expr=blk.pump.efficiency_fluid
        == blk.pump.efficiency_eq_cubed * flow**3
        + blk.pump.efficiency_eq_squared * flow**2
        + blk.pump.efficiency_eq_linear * flow
        + blk.pump.efficiency_eq_constant,
        doc="Efficiency surrogate equation",
    )
    blk.pump.efficiency_pump.bounds = (0, 1)

    if stage_num == 1:
        efficiency_motor = 0.962
        efficiency_vfd = 0.97
    elif stage_num == 2:
        efficiency_motor = (
            0.938  # Likely these numbers are lower (due to elevated temperatures ?)
        )
        efficiency_vfd = 0.97
    elif stage_num == 3:
        # Making same as stage 2 for now
        efficiency_motor = 0.938
        efficiency_vfd = 0.97

    blk.pump.efficiency_motor = Param(
        initialize=efficiency_motor,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Efficiency of motor and VFD",
    )

    blk.pump.efficiency_vfd = Param(
        initialize=efficiency_vfd,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Efficiency of motor and VFD",
    )

    blk.pump.efficiency_electrical = Constraint(
        expr=blk.pump.efficiency_pump[0]
        == blk.pump.efficiency_fluid
        * blk.pump.efficiency_motor
        * blk.pump.efficiency_vfd
    )

    # Add Arcs
    blk.feed_to_unit = Arc(source=blk.feed.outlet, destination=blk.pump.inlet)
    blk.unit_to_product = Arc(source=blk.pump.outlet, destination=blk.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_pump_op_conditions(blk):
    Pout = get_config_value(
        blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{blk.stage_num}"
    )
    blk.pump.control_volume.properties_out[0].pressure.fix(Pout)


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
    blk.pump.control_volume.properties_in[0].flow_vol  # Touching

    # Scaling properties
    m = blk.model()
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def add_pump_scaling(blk):
    # Properties
    set_scaling_factor(blk.pump.work_mechanical[0], 1e-3)
    set_scaling_factor(blk.pump.efficiency_pump, 1e1)


def initialize_pump(blk):
    # Touch Properties that are needed for later
    blk.feed.properties[0].flow_vol_phase["Liq"]
    blk.product.properties[0].flow_vol_phase["Liq"]

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.pump.initialize()
    propagate_state(blk.unit_to_product)
    blk.product.initialize()


def report_pump(blk, w=30):
    title = "Pump Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    total_flow = blk.pump.control_volume.properties_in[0].flow_vol_phase["Liq"]
    deltaP = blk.pump.deltaP[0]
    work = blk.pump.work_mechanical[0]
    print(
        f'{f"Total Flow Rate (MGD)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.Mgallons /pyunits.day)):<{w}.3f}{"MGD"}'
    )
    print(f'{f"Total Flow Rate (m3/s)":<{w}s}{value(total_flow):<{w}.3e}{"m3/s"}')
    print(
        f'{f"Total Flow Rate (gpm)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Pressure Change (psi)":<{w}s}{value(pyunits.convert(deltaP, to_units=pyunits.psi)):<{w}.3e}{"psi"}'
    )
    print(
        f'{f"Head (ft)":<{w}s}{value(pyunits.convert(deltaP, to_units=pyunits.ftH2O)):<{w}.3e}{"ft"}'
    )
    print(
        f'{f"Pressure Change (bar)":<{w}s}{value(pyunits.convert(deltaP,to_units=pyunits.bar)):<{w}.3e}{"bar"}'
    )
    print(f'{f"Work Mech. (J)":<{w}s}{value(work):<{w}.3e}{"Joules"}')
    print(
        f'{f"Work Mech. (kW)":<{w}s}{value(pyunits.convert(work, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    print(f'{f"Efficiency (-)":<{w}s}{value(blk.pump.efficiency_pump[0]):<{w}.3f}{"-"}')


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
    stage_num = 1
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
