import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    TransformationFactory,
    Objective,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
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
    Product,
    Feed,
    StateJunction,
    Mixer,
    Separator,
    EnergySplittingType,
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

from watertap.flowsheets.flex_desal.wrd.components.ro_system import (
    load_config,
    get_config_value,
)

from watertap.core import Database


# TODO:
# 1. Unfix the variable energy_electric_flow_vol_inlet


def build_UF_pumps_system(split_fractions, config=None):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = NaClParameterBlock()

    m.fs.UF_pumps = FlowsheetBlock(dynamic=False)
    build_UF_pumps(
        m.fs.UF_pumps, m.fs.properties, split_fractions=split_fractions, config=config
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_UF_pumps(
    blk, prop_package, number_trains=3, split_fractions=None, config=None
) -> None:

    print(f'\n{"=======> BUILDING ULTRAFILTRATION SYSTEM <=======":^60}\n')

    current_script_path = os.path.abspath(__file__)
    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)
    # Get the parent directory of the current directory (one folder prior)
    parent_directory = os.path.dirname(current_directory)

    print(f"Parent directory: {parent_directory}")

    if config is None:
        config = parent_directory + "/meta_data/wrd_uf_pumps_inputs.yaml"
    else:
        config = "/Users/mhardika/Documents/watertap/watertap/watertap/flowsheets/flex_desal/wrd/meta_data/wrd_uf_pumps_inputs.yaml"
    blk.config_data = load_config(config)

    blk.number_trains = number_trains

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.UF_properties

    blk.feed_in = StateJunction(property_package=prop_package)
    blk.feed_out = StateJunction(property_package=prop_package)

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
        source=blk.feed_in.outlet,
        destination=blk.feed_splitter.inlet,
    )

    blk.pump_mixer_to_feed_out = Arc(
        source=blk.pump_outlet_mixer.outlet,
        destination=blk.feed_out.inlet,
    )


def set_UF_pumps_inlet_conditions(blk, Qin=0.618, Cin=0.542):
    """
    Set the operation conditions for the UF pumps
    """
    Qin = (Qin) * pyunits.m**3 / pyunits.s  # Feed flow rate in m3/s
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed_in.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        feed_mass_flow_water
    )
    blk.feed_in.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        feed_mass_flow_salt
    )
    blk.feed_in.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed_in.properties[0].pressure.fix(101325 * pyunits.Pa)  # 1 bar


def set_UF_pump_op_conditions(blk):
    # Set pump operating conditions
    for i in range(1, blk.number_trains + 1):
        pump = blk.find_component(f"pump_{i}")

        pump.control_volume.properties_out[0].pressure.fix(
            get_config_value(
                blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{i}"
            )
        )

        pump.efficiency_pump.fix(
            get_config_value(blk.config_data, "pump_efficiency", "pumps", f"pump_{i}")
        )


def init_UF_pumps(blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    print(
        "\n\n-------------------- INITIALIZING ULTRAFILTRATION --------------------\n\n"
    )

    print(f"UF Pumps Degrees of Freedom: {degrees_of_freedom(blk)}")
    print("\n\n")

    blk.feed_in.initialize(optarg=optarg)

    propagate_state(blk.feed_to_feed_splitter)

    blk.feed_splitter.initialize()

    for i in range(1, blk.number_trains + 1):

        splitter_out = blk.feed_splitter.find_component(f"pump_{i}_feed")

        # Propagate state to each pump
        pump = blk.find_component(f"pump_{i}")
        splitter_to_pump = blk.find_component(f"splitter_to_pump_{i}_connect")
        propagate_state(splitter_to_pump)

        pump.initialize()

        propagate_state(blk.find_component(f"pump_{i}_to_mixer_connect"))
        blk.pump_outlet_mixer.initialize()

        propagate_state(blk.pump_mixer_to_feed_out)
        blk.feed_out.initialize(optarg=optarg)


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


if __name__ == "__main__":

    split_fractions = [
        0.4,
        0.4,
        0.2,
    ]  # Based on ratio of pump capacity to total capacity

    m = build_UF_pumps_system(split_fractions=split_fractions)

    set_UF_pumps_inlet_conditions(m.fs.UF_pumps)
    set_UF_pump_op_conditions(m.fs.UF_pumps)

    init_UF_pumps(m.fs.UF_pumps)
    solve(m)

    add_costing(m)
    add_UF_pumps_costing(m, m.fs.UF_pumps)
    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    solve(m)

    print("Degrees of Freedom: ", degrees_of_freedom(m))

    # Electricity consumption
    m.fs.costing.aggregate_flow_electricity.display()

    # Print work mechanical of each pump
    for i in range(1, (m.fs.UF_pumps.number_trains + 1)):
        pump = m.fs.UF_pumps.find_component(f"pump_{i}")
        electricity = pyunits.convert(pump.work_mechanical[0], to_units=pyunits.kW)
        print(
            f"Pump {i} Work Mechanical: {value(electricity) :.2f} {pyunits.get_units(electricity)}"
        )
