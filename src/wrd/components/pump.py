from pyomo.environ import (
    ConcreteModel,
    Var,
    Param,
    Constraint,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
    Block,
)
from pyomo.network import Arc
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import StateJunction, Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor

from watertap.costing import WaterTAPCosting
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.unit_models.pressure_changer import Pump
from watertap.core.solvers import get_solver

from wrd.utilities import (
    load_config,
    get_config_value,
    get_config_file,
    ft_head_to_psi,
    psi_to_ft_head,
)
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


def build_system(
    stage_num=1, uf=False, Qin=None, head=None, file="wrd_inputs_8_19_21.yaml"
):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)
    m.fs.pump = FlowsheetBlock(dynamic=False)
    build_pump(
        m.fs.pump,
        stage_num=stage_num,
        uf=uf,
        file=file,
        prop_package=m.fs.properties,
        Qin=Qin,
        head=head,
    )
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


def find_pump_speed(blk, stage_num=1, uf=False, Qin=None, head=None):
    # Create variable for pump speed and reference (100% speed) flow and head

    blk.unit.eff.speed = Var(
        initialize=1,
        units=pyunits.dimensionless,
        bounds=(0, 1.02),
        doc="Pump speed ratio (actual speed / maximum speed)",
    )

    blk.unit.eff.ref_head = Var(
        initialize=1,
        units=pyunits.feet,
        bounds=(0, None),  # Should actually have bounds
        doc="Pump reference head at 100% speed (ft)",
    )

    blk.unit.eff.ref_flow = Var(
        initialize=0.1,
        units=pyunits.m**3 / pyunits.s,
        bounds=(0, None),  # Should actually have bounds
        doc="Pump reference flow at 100% speed (m3/s)",
    )

    if head is None:
        if uf:
            head = psi_to_ft_head(
                get_config_value(
                    blk.config_data, "pump_outlet_pressure", "uf_pumps", "pump"
                )
                - 14.5 * pyunits.psi
            )  # Assuming atmopheric suct. pressure
        else:
            head = psi_to_ft_head(
                get_config_value(
                    blk.config_data,
                    "pump_outlet_pressure",
                    "ro_pumps",
                    f"pump_stage_{stage_num}",
                )
                - get_config_value(
                    blk.config_data,
                    "pump_suction_pressure",
                    "ro_pumps",
                    f"pump_stage_{stage_num}",
                )
            )

    flow = pyunits.convert(Qin, to_units=pyunits.m**3 / pyunits.s)

    blk.unit.eff.head = Param(
        initialize=head,
        units=pyunits.feet,
        mutable=True,
        doc="Pump pressure differential as head (ft)",
    )

    # Possibly redundant
    blk.unit.eff.flow = Param(
        initialize=flow,
        units=pyunits.m**3 / pyunits.s,
        mutable=True,
        doc="Feed Flowrate",
    )

    # EQ 1
    # the pump affinity laws into Constraints
    blk.unit.eff.eq_head_affinity_law = Constraint(
        expr=blk.unit.eff.head == blk.unit.eff.ref_head * blk.unit.eff.speed**2,
        doc="Pump head affinity law equation",
    )
    # EQ 2
    blk.unit.eff.eq_flow_affinity_law = Constraint(
        expr=blk.unit.eff.flow == blk.unit.eff.ref_flow * blk.unit.eff.speed,
        doc="Pump flow affinity law equation",
    )
    # Load head curve for 100% speed
    if uf:
        b_0 = 323.88
        b_1 = -403.68
        b_2 = 1449.8
        b_3 = -6297.6
    elif stage_num == 1:
        b_0 = 374.73
        b_1 = -1347.1
        b_2 = 8953.8
        b_3 = -26539
    elif stage_num == 2:
        b_0 = 100.72
        b_1 = -189.57
        b_2 = 4535.7
        b_3 = -87520
    else:
        # Still missing TSRO pump curves
        b_0 = 0
        b_1 = 0
        b_2 = 0
        b_3 = 0

    # Create parameters for the fit for 100% speed
    blk.unit.eff.ref_head_constant = Param(
        initialize=b_0,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Constant term of Efficiency equation",
    )

    blk.unit.eff.ref_head_linear_coeff = Param(
        initialize=b_1,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -1,
        doc="Linear term of Efficiency equation",
    )

    blk.unit.eff.ref_head_squared_coeff = Param(
        initialize=b_2,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -2,
        doc="Squared term of Efficiency equation",
    )

    blk.unit.eff.ref_head_cubed_coeff = Param(
        initialize=b_3,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -3,
    )

    # EQ 3
    blk.unit.eff.ref_head_surr = Constraint(
        expr=blk.unit.eff.ref_head
        == blk.unit.eff.ref_head_cubed_coeff * blk.unit.eff.ref_flow**3
        + blk.unit.eff.ref_head_squared_coeff * blk.unit.eff.ref_flow**2
        + blk.unit.eff.ref_head_linear_coeff * blk.unit.eff.ref_flow
        + blk.unit.eff.ref_head_constant,
        doc="Head surrogate equation",
    )

    # calculate_variable_from_constraint(blk.unit.eff.speed, blk.unit.eff.eq_head_affinity_law) # Is that the right eq?
    # solver.solve(blk.unit.eff)
    # print(f"Calculated pump speed for stage {stage_num}: {value(blk.unit.eff.speed)}")


def set_pump_efficiency(blk, stage_num=1, uf=False, Qin=None, head=None):

    # Creating a subblock for all the efficiency related vars, param, and constraints. That way, they can be solved without solve whole pump for trouble shooting. Can remove if not useful later.
    blk.unit.eff = Block()
    blk.unit.eff.efficiency_fluid = Var(
        initialize=1,
        units=pyunits.dimensionless,
        bounds=(0, 1),
        doc="Pump efficiency from pump curves",
    )
    # Load Values for efficiency surrogate model
    if uf:
        # Below are estimated for 100% speed
        a_0 = 0.0677
        a_1 = 5.357
        a_2 = -4.475
        a_3 = -19.578
    elif stage_num == 1:
        a_0 = 0.389
        a_1 = -0.535
        a_2 = 41.373
        a_3 = -138.820
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

    # Create Variables for max speed efficiency curve
    blk.unit.eff.efficiency_constant = Param(
        initialize=a_0,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Constant term of Efficiency equation",
    )

    blk.unit.eff.efficiency_linear_coeff = Param(
        initialize=a_1,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -1,
        doc="Linear term of Efficiency equation",
    )

    blk.unit.eff.efficiency_squared_coeff = Param(
        initialize=a_2,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -2,
        doc="Squared term of Efficiency equation",
    )

    blk.unit.eff.efficiency_cubed_coeff = Param(
        initialize=a_3,
        mutable=True,
        units=(pyunits.m**3 / pyunits.s) ** -3,
        doc="Cubed term of Efficiency equation",
    )

    find_pump_speed(blk, stage_num=stage_num, uf=uf, Qin=Qin, head=head)

    # ref_flow = flow at 100% speed with the same efficiency
    # ref_flow = blk.feed.properties[0].flow_vol_phase["Liq"] / blk.unit.eff.speed

    blk.unit.eff.eq_efficiency_surr = Constraint(
        expr=blk.unit.eff.efficiency_fluid
        == blk.unit.eff.efficiency_cubed_coeff * blk.unit.eff.ref_flow**3
        + blk.unit.eff.efficiency_squared_coeff * blk.unit.eff.ref_flow**2
        + blk.unit.eff.efficiency_linear_coeff * blk.unit.eff.ref_flow
        + blk.unit.eff.efficiency_constant,
        doc="Efficiency surrogate equation",
    )
    blk.unit.efficiency_pump.bounds = (0, 1)  # Is this needed?
    assert degrees_of_freedom(blk.unit.eff) == 0
    solver.solve(blk.unit.eff)
    # print(f"Calculated pump speed for stage {stage_num}: {value(blk.unit.eff.speed)}")


def build_pump(
    blk,
    stage_num=1,
    file="wrd_inputs_8_19_21.yaml",
    prop_package=None,
    uf=False,
    Qin=None,
    head=None,
):

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.ro_properties

    blk.config_data = load_config(get_config_file(file))
    blk.stage_num = stage_num

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = Pump(property_package=prop_package)

    blk.product = StateJunction(property_package=prop_package)
    set_pump_efficiency(blk, stage_num=stage_num, uf=uf, Qin=Qin, head=head)

    # Create variable for the efficiency from the pump curves
    blk.unit.efficiency_motor = Param(
        initialize=0.95,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Efficiency of motor and VFD",
    )

    blk.unit.efficiency_vfd = Param(
        initialize=0.97,
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
            * blk.unit.eff.efficiency_fluid
        )
        - blk.unit.efficiency_loss
    )

    # Add Arcs
    blk.feed_to_unit = Arc(source=blk.feed.outlet, destination=blk.unit.inlet)
    blk.unit_to_product = Arc(source=blk.unit.outlet, destination=blk.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_pump_op_conditions(blk, uf=False, head=None, Pin=14.5):

    if head is None:
        if uf:
            # All the pumps are assumed to have the same outlet pressure for UF pumps because they collect in a header
            Pout = get_config_value(
                blk.config_data, "pump_outlet_pressure", "uf_pumps", f"pump"
            )
        else:
            Pout = get_config_value(
                blk.config_data,
                "pump_outlet_pressure",
                "ro_pumps",
                f"pump_stage_{blk.stage_num}",
            )
            print(
                f"Setting pump {blk.stage_num} operating conditions, Pout = {value(Pout)} psi"
            )
    else:
        head = head * pyunits.feet
        Pout = ft_head_to_psi(head) + Pin * pyunits.psi

    blk.unit.control_volume.properties_out[0].pressure.fix(Pout)


def set_inlet_conditions(m, Qin=None, Cin=0.5, Tin=302, Pin=14.5):
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): Pin * pyunits.psi,
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


def add_pump_costing(blk, costing_package=None):

    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_package)
    # Only want to cost opex (electricity)
    costing_package.high_pressure_pump.cost.fix(0)


def report_pump(blk, w=30, add_costing=False):
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
    print(
        f'{f"Pump Speed Ratio (%)":<{w}s}{100*value(blk.unit.eff.speed):<{w}.3f}{"%"}'
    )

    print(f'{f"Efficiency (-)":<{w}s}{value(blk.unit.efficiency_pump[0]):<{w}.3f}{"-"}')
    if add_costing:
        m = blk.model()
        # Is SEC not appearing on m.fs.costing.display a known issue?
        SEC = m.fs.costing.SEC
        print(
            f'{f"Specific Energy (SEC)":<{w}s}{value(pyunits.convert(SEC, to_units=pyunits.kWh / pyunits.m**3)):<{w}.3f}{"kWh/m3"}'
        )


def main(
    Qin=None,
    head=None,  # Entering head value will override the Pout value in the yaml
    Cin=0.5,
    Tin=302,
    Pin=14.5,  # psi
    stage_num=1,
    uf=False,
    file="wrd_inputs_8_19_21.yaml",
    add_costing=True,
):
    # Handling loading Qin here b/c it's an input to both build_pump and set_inlet_conditions
    if Qin is None:
        if uf:
            Qin = pyunits.convert(
                get_config_value(blk.config_data, "pump_flowrate", "uf_pumps", "pump"),
                to_units=pyunits.gal / pyunits.minute,
            )
        else:
            # load the flowrate and pressure head  WATCH OUT FOR UNITS. These are knowns need to calc above values.
            Qin = pyunits.convert(
                get_config_value(
                    blk.config_data,
                    "pump_flowrate",
                    "ro_pumps",
                    f"pump_stage_{stage_num}",
                ),
                to_units=pyunits.gal / pyunits.minute,
            )
    else:
        Qin = Qin * pyunits.gal / pyunits.minute

    m = build_system(stage_num=stage_num, uf=uf, file=file, Qin=Qin, head=head)
    add_pump_scaling(m.fs.pump)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin, Tin=Tin, Pin=Pin)
    set_pump_op_conditions(m.fs.pump, head=head, Pin=Pin, uf=uf)

    if add_costing:
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
    report_pump(m.fs.pump, add_costing=add_costing)

    return m


if __name__ == "__main__":

    # August 19, 2021 Data
    # Stage 1
    # m = main()
    # Testing at a lower speed
    m = main(
        Qin=3894,
        head=None,
        Cin=1.2,
        Tin=302,
        Pin=14.5,
        stage_num=None,
        file="wrd_inputs_8_19_21.yaml",
        uf=True,
    )
    # Stage 2
    # m = main(Qin=1029, Pin=131.2 * pyunits.psi, stage_num=2)
    # Stage 3
    # m = main(Qin=384, Pin=(112.6 - 41.9) * pyunits.psi, stage_num=3)
