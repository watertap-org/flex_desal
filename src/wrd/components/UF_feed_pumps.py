from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    TransformationFactory,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
from idaes.core import (
    MomentumBalanceType,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import (
    StateJunction,
    Mixer,
    Separator,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import *
from watertap.core.util.initialization import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.unit_models.zero_order.ultra_filtration_zo import UltraFiltrationZO
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.unit_models.pressure_changer import Pump

from wrd.utilities import (
    get_config_file,
    load_config,
    get_config_value,
)

# TODO:
# 1. Unfix the variable energy_electric_flow_vol_inlet


def build_system(split_fractions):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.UF_pumps = FlowsheetBlock(dynamic=False)
    build_UF_pumps(
        m.fs.UF_pumps,
        m.fs.properties,
        split_fractions=split_fractions,
    )
    return m


def build_UF_pumps(blk, prop_package, split_fractions=None, date="8_19_21") -> None:

    # print(f'\n{"=======> BUILDING ULTRAFILTRATION SYSTEM <=======":^60}\n')

    config_file_name = get_config_file("wrd_uf_pumps_inputs_" + date + ".yaml")
    blk.config_data = load_config(config_file_name)
    number_trains = get_config_value(blk.config_data, "number_trains", "pumps")
    blk.number_trains = number_trains  # Could be moved to config / yaml
    assert (
        len(split_fractions) == number_trains
    ), "Length of split fractions must equal number of trains"
    assert sum(split_fractions[0:number_trains]) == 1

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.ro_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)

    # Splits the feed to multiple UF feed pumps
    blk.feed_splitter = Separator(
        property_package=prop_package,
        outlet_list=[f"pump_{i+1}_feed" for i in range(number_trains)],
    )

    # Combines the outlet of all UF feed pumps
    blk.pump_outlet_mixer = Mixer(
        property_package=prop_package,
        inlet_list=[f"pump_{i+1}_to_mixer" for i in range(number_trains)],
        energy_mixing_type=MixingType.extensive,
        momentum_mixing_type=MomentumMixingType.minimize,
    )

    if split_fractions is None:
        split_fractions = [1 / number_trains for i in range(number_trains)]

    for i in range(1, (blk.number_trains + 1)):
        # Add pump for each train
        blk.add_component(
            f"pump_{i}",
            Pump(property_package=prop_package),
        )

        # Calculate flow to each pump
        blk.feed_splitter.split_fraction[0, f"pump_{i}_feed"].set_value(
            split_fractions[i - 1]
        )
        if i != 1:
            blk.feed_splitter.split_fraction[0, f"pump_{i}_feed"].fix()

        # Connect splitter outlet to pump inlet
        splitter_outlet = blk.feed_splitter.find_component(f"pump_{i}_feed")
        blk.add_component(
            f"splitter_to_pump_{i}_connect",
            Arc(
                source=splitter_outlet,
                destination=blk.find_component(f"pump_{i}").inlet,
            ),
        )

        # Connect pump outlet to mixer inlet
        blk.add_component(
            f"pump_{i}_to_mixer_connect",
            Arc(
                source=blk.find_component(f"pump_{i}").outlet,
                destination=blk.pump_outlet_mixer.find_component(f"pump_{i}_to_mixer"),
            ),
        )

    blk.feed_to_feed_splitter = Arc(
        source=blk.feed.outlet,
        destination=blk.feed_splitter.inlet,
    )

    blk.pump_mixer_to_product = Arc(
        source=blk.pump_outlet_mixer.outlet,
        destination=blk.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_inlet_conditions(blk, Qin=None, Cin=None):
    """
    Set the operation conditions for the UF pumps.

    If Qin or Cin are not provided they are read from the block's
    `config_data` (using the same keys/structure as the pump inlet helper).
    This mirrors `set_inlet_conditions` in `wrd/components/pump.py`.
    """
    # Read from config if not provided. Config values are expected to
    # already carry units consistent with the rest of the model (as in
    # `set_inlet_conditions` in `pump.py`). Use pump_1 as the default key.
    if Qin is None:
        Qin = get_config_value(blk.config_data, "feed_flow_water", "feed_stream")

    if Cin is None:
        Cin = get_config_value(
            blk.config_data, "feed_conductivity", "feed_stream"
        ) * get_config_value(
            blk.config_data, "feed_conductivity_conversion", "feed_stream"
        )

    Pin = get_config_value(blk.config_data, "feed_pressure", "feed_stream")

    # Approximate density of water
    rho = 1000 * pyunits.kg / pyunits.m**3
    # Calculate mass flows from volumetric flow and concentration
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    blk.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed.properties[0].pressure.fix(Pin)
    # Touching volumetric flow variables for initialization
    blk.feed.properties[0].flow_vol_phase["Liq"]

    # Scaling defaults on the top-level property block (match pump behavior)
    m = blk.model()
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def set_UF_pumps_op_conditions(blk):
    # Set pump operating conditions
    for i in range(1, blk.number_trains + 1):
        pump = blk.find_component(f"pump_{i}")

        pump.control_volume.properties_out[0].pressure.fix(
            get_config_value(
                blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{i}"
            )
        )
        pump.control_volume.properties_in[0].flow_vol_phase[
            "Liq"
        ]  # Touching so it can be used in later EQ
        # Create variable for the efficiency from the pump curves on each pump
        pump.efficiency_fluid = Var(
            initialize=0.5,
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Efficiency from pump curves",
        )

        # Load Values for surrogate model for 100% pump speed.
        a_0 = 0.0677
        a_1 = 5.357
        a_2 = -4.475
        a_3 = -19.578

        # Create Params for simple "surrogate"
        pump.efficiency_eq_constant = Param(
            initialize=a_0,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Constant term of Efficiency equation",
        )

        pump.efficiency_eq_linear = Param(
            initialize=a_1,
            mutable=True,
            units=(pyunits.m**3 / pyunits.s) ** -1,
            doc="Linear term of Efficiency equation",
        )

        pump.efficiency_eq_squared = Param(
            initialize=a_2,
            mutable=True,
            units=(pyunits.m**3 / pyunits.s) ** -2,
            doc="Squared term of Efficiency equation",
        )

        pump.efficiency_eq_cubed = Param(
            initialize=a_3,
            mutable=True,
            units=(pyunits.m**3 / pyunits.s) ** -3,
            doc="Cubed term of Efficiency equation",
        )

        flow = pump.control_volume.properties_in[0].flow_vol_phase["Liq"]

        pump.efficiency_surr_eq = Constraint(
            expr=pump.efficiency_fluid
            == pump.efficiency_eq_cubed * flow**3
            + pump.efficiency_eq_squared * flow**2
            + pump.efficiency_eq_linear * flow
            + pump.efficiency_eq_constant,
            doc="Efficiency surrogate equation",
        )

        # Ensure pump efficiency bounds/value are sensible
        pump.efficiency_pump.bounds = (0, 1)

        pump.efficiency_motor = Param(
            initialize=0.95,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Efficiency of motor",
        )

        pump.efficiency_vfd = Param(
            initialize=0.95,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Efficiency of VFD",
        )

        pump.efficiency_electrical = Constraint(
            expr=pump.efficiency_pump[0]
            == pump.efficiency_motor * pump.efficiency_vfd * pump.efficiency_fluid
        )


def add_UF_pump_scaling(blk):
    for i in range(1, blk.number_trains + 1):
        pump = blk.find_component(f"pump_{i}")
        set_scaling_factor(pump.work_mechanical[0], 1e-3)


def init_UF_pumps(blk, verbose=True, solver=None):

    if solver is None:
        solver = get_solver()
    # Why the optarg here? This approach is not used in any other component
    optarg = solver.options
    # print(
    #     "\n\n-------------------- INITIALIZING ULTRAFILTRATION --------------------\n\n"
    # )

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_feed_splitter)

    blk.feed_splitter.initialize()

    for i in range(1, blk.number_trains + 1):

        # Propagate state to each pump
        pump = blk.find_component(f"pump_{i}")
        splitter_to_pump = blk.find_component(f"splitter_to_pump_{i}_connect")
        propagate_state(splitter_to_pump)
        pump.initialize()
        propagate_state(blk.find_component(f"pump_{i}_to_mixer_connect"))
        blk.pump_outlet_mixer.initialize()
        propagate_state(blk.pump_mixer_to_product)

    blk.product.initialize(optarg=optarg)


def add_UF_pumps_costing(m, blk, costing_blk=None):
    if costing_blk is None:
        costing_blk = m.fs.costing

    for i in range(1, (blk.number_trains + 1)):
        pump = blk.find_component(f"pump_{i}")
        pump.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_blk)


def add_costing(m):
    m.fs.costing = ZeroOrderCosting()


def solve(model, solver=None, tee=True, raise_on_failure=True):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results


def print_UF_costing_breakdown(blk, debug=False):
    print(f"\n\n-------------------- UF Costing Breakdown --------------------\n")
    print(f'{"UF Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}')

    if debug:
        print(blk.unit.costing.display())


def report_UF_pumps(blk, w=30):
    title = "UF Feed Pumps Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    # Print work mechanical of each pump
    for i in range(1, (blk.number_trains + 1)):
        pump = blk.find_component(f"pump_{i}")
        electricity = pyunits.convert(pump.work_mechanical[0], to_units=pyunits.kW)
        flow = pump.control_volume.properties_in[0].flow_vol_phase["Liq"]
        deltaP = pump.deltaP[0]
        print("\n")
        print(
            f'{f"Pressure Change (psi)":<{w}s}{value(pyunits.convert(deltaP,to_units=pyunits.psi)):<{w}.3e}{"psi"}'
        )
        print(
            f'{f"Pump {i} Flow Rate (gpm)":<{w}s}{value(pyunits.convert(flow, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(f'{f"Efficiency (-)":<{w}s}{value(pump.efficiency_pump[0]):<{w}.3f}{"-"}')
        print(
            f'{f"Pump {i} Work Mech. (kW)":<{w}s}{value(pyunits.convert(electricity, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
        )


if __name__ == "__main__":
    # Probably will want to move split fraction into the yaml
    split_fractions = [1]# Based on ratio of pump capacity to total capacity

    m = build_system(split_fractions=split_fractions)
    assert_units_consistent(m)
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_inlet_conditions(m.fs.UF_pumps)
    set_UF_pumps_op_conditions(m.fs.UF_pumps)
    print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    add_UF_pump_scaling(m.fs.UF_pumps)
    calculate_scaling_factors(m)
    init_UF_pumps(m.fs.UF_pumps)
    solve(m)
    report_UF_pumps(m.fs.UF_pumps)

    # Not doing any costing for the moment
    # add_costing(m)
    # add_UF_pumps_costing(m, m.fs.UF_pumps)
    # m.fs.costing.cost_process()
    # m.fs.costing.initialize()
    # solve(m)

    # Electricity consumption
    # m.fs.costing.aggregate_flow_electricity.display()
