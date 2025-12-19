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

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import StateJunction, Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor

from watertap.costing import WaterTAPCosting
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
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
    "add_pump_costing",
]

solver = get_solver()


def build_system(stage_num=1, file="wrd_inputs_8_19_21.yaml"):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)
    m.fs.pump = FlowsheetBlock(dynamic=False)
    build_pump(m.fs.pump, stage_num=stage_num, file=file, prop_package=m.fs.properties)

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)

    # Arcs to connect the unit models
    m.fs.feed_to_pump = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.pump.feed.inlet,
    )
    m.fs.pump_to_product = Arc(
        source=m.fs.pump.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def build_pump(blk, stage_num=1, file="wrd_inputs_8_19_21.yaml", prop_package=None,uf=False):

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.ro_properties

    blk.config_data = load_config(get_config_file(file))
    blk.stage_num = stage_num

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = Pump(property_package=prop_package)

    blk.product = StateJunction(property_package=prop_package)

    # Create variable for the efficiency from the pump curves
    blk.unit.efficiency_fluid = Var(
        initialize=0.7,
        units=pyunits.dimensionless,
        bounds=(0, 1),
        doc="Efficiency from pump curves",
    )

    # Load Values for surrogate model
    if uf:
        a_0 = 0.0677
        a_1 = 5.357
        a_2 = -4.475
        a_3 = -19.578
    elif stage_num == 1:
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
        # Still missing TSRO pump curves
        a_0 = 0.067
        a_1 = 21.112
        a_2 = -133.157
        a_3 = -234.386

    # Create Variables for simple "surrogate"
    blk.unit.efficiency_constant = Param(
        initialize=a_0,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Constant term of Efficiency equation",
    )

    blk.unit.efficiency_linear_coeff = Param(
        initialize=a_1,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -1,
        doc="Linear term of Efficiency equation",
    )

    blk.unit.efficiency_squared_coeff = Param(
        initialize=a_2,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -2,
        doc="Squared term of Efficiency equation",
    )

    blk.unit.efficiency_cubed_coeff = Param(
        initialize=a_3,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -3,
        doc="Cubed term of Efficiency equation",
    )

    flow = blk.feed.properties[0].flow_vol_phase["Liq"]

    blk.unit.eq_efficiency_surr = Constraint(
        expr=blk.unit.efficiency_fluid
        == blk.unit.efficiency_cubed_coeff * flow**3
        + blk.unit.efficiency_squared_coeff * flow**2
        + blk.unit.efficiency_linear_coeff * flow
        + blk.unit.efficiency_constant,
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

    blk.unit.efficiency_loss = Param(
        initialize=0,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Loss factor due to heat, age, wear, etc.",
    )

    blk.unit.eq_efficiency_electrical = Constraint(
        expr=blk.unit.efficiency_pump[0]
        == (
            blk.unit.efficiency_motor
            * blk.unit.efficiency_vfd
            * blk.unit.efficiency_fluid
        )
        - blk.unit.efficiency_loss
    )

    # Add Arcs
    blk.feed_to_unit = Arc(source=blk.feed.outlet, destination=blk.unit.inlet)
    blk.unit_to_product = Arc(source=blk.unit.outlet, destination=blk.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_pump_op_conditions(blk, uf=False):
    if uf:
        Pout = get_config_value(
            blk.config_data, "pump_outlet_pressure", "uf_pumps", f"pump_{blk.stage_num}"
        )
    else:
        Pout = get_config_value(
            blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{blk.stage_num}"
        )
        print(
            f"Setting pump {blk.stage_num} operating conditions, Pout = {value(Pout)} psi"
        )
    blk.unit.control_volume.properties_out[0].pressure.fix(Pout)


def set_inlet_conditions(m, Qin=2637, Cin=0.5, Tin=302, Pin=101325):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): Pin,
            ("temperature", None): Tin,
        },
        hold_state=True,
    )


def add_pump_scaling(blk):
    set_scaling_factor(blk.unit.work_mechanical[0], 1e-3)


def initialize_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_pump)

    initialize_pump(m.fs.pump)

    propagate_state(m.fs.pump_to_product)
    m.fs.product.initialize()


def initialize_pump(blk):

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
    work = blk.unit.work_mechanical[0]
    pin = blk.unit.control_volume.properties_in[0].pressure
    deltaP = blk.unit.deltaP[0]
    pout = blk.unit.control_volume.properties_out[0].pressure
    print(
        f'{f"Inlet Flow":<{w}s}{value(pyunits.convert(flow_in, to_units=pyunits.gallons /pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    # print(f'{f"∆P (Pa)":<{w}s}{value(deltaP):<{w}.3e}{"Pa"}')
    print(
        f'{f"Inlet Pressure":<{w}s}{value(pyunits.convert(pin, to_units=pyunits.psi)):<{w}.3f}{"psi"}'
    )
    print(
        f'{f"∆P":<{w}s}{value(pyunits.convert(deltaP, to_units=pyunits.psi)):<{w}.3f}{"psi"}'
    )
    print(
        f'{f"Outlet Pressure":<{w}s}{value(pyunits.convert(pout, to_units=pyunits.psi)):<{w}.3f}{"psi"}'
    )
    print(
        f'{f"Work Mech. (kW)":<{w}s}{value(pyunits.convert(work, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    print(f'{f"Efficiency (-)":<{w}s}{value(blk.unit.efficiency_pump[0]):<{w}.3f}{"-"}')


def add_pump_costing(blk, costing_package=None):

    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    # blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_package)
    costing_package.cost_flow(
        pyunits.convert(blk.unit.work_mechanical[0], to_units=pyunits.kW), "electricity"
    )


def main(
    Qin=2637,
    Cin=0.5,
    Tin=302,
    Pin=101325,
    stage_num=1,
    file="wrd_inputs_8_19_21.yaml",
):

    m = build_system(stage_num=stage_num, file=file)
    add_pump_scaling(m.fs.pump)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin, Tin=Tin, Pin=Pin)
    set_pump_op_conditions(m.fs.pump)

    add_pump_costing(m.fs.pump)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"],
        name="SEC",
    )

    initialize_system(m)
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_pump(m.fs.pump)

    return m


if __name__ == "__main__":

    # August 19, 2021 Data
    # Stage 1
    m = main()
    # Stage 2
    m = main(Qin=1029, Pin=131.2 * pyunits.psi, stage_num=2)
    # Stage 3
    m = main(Qin=384, Pin=(112.6 - 41.9) * pyunits.psi, stage_num=3)
    m.fs.costing.SEC.display()
