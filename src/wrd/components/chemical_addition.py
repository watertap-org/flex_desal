import pathlib
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent
from pyomo.util.calc_var_value import calculate_variable_from_constraint as cvc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.core.solvers import get_solver
from watertap.core import Database
from watertap.unit_models.zero_order import ChemicalAdditionZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.util.initialization import *


def build_system():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.db = Database(dbpath="src/wrd/meta_data")
    m.fs.properties = WaterParameterBlock(solute_list=["tds", "tss"])

    m.fs.chem_addition = FlowsheetBlock(dynamic=False)

    build_chem_addition(m.fs.chem_addition, "ammonia", m.fs.properties)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_chem_addition(blk, chemical_name, prop_package=None) -> None:

    print(f'\n{"=======> BUILDING CHEMICAL ADDITION SYSTEM <=======":^60}\n')

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)

    blk.unit = ChemicalAdditionZO(
        property_package=prop_package,
        database=m.db,
        process_subtype=chemical_name,
    )
    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )


def set_system_conditions(blk):
    blk.feed.properties[0.0].flow_mass_comp["H2O"].fix(171.37)
    blk.feed.properties[0.0].flow_mass_comp["tds"].fix(600)
    blk.feed.properties[0.0].flow_mass_comp["tss"].fix(5.22e-6)


def set_chem_addition_scaling(blk, calc_blk_scaling_factors=False):

    set_scaling_factor(blk.unit.chemical_dosage, 0.1)
    set_scaling_factor(blk.unit.solution_density, 1e-3)
    set_scaling_factor(blk.unit.chemical_flow_vol, 1e6)
    set_scaling_factor(blk.unit.electricity, 1e4)

    # Calculate scaling factors only for chem addition block if in full case study flowsheet
    # so we don't prematurely set scaling factors
    if calc_blk_scaling_factors:
        calculate_scaling_factors(blk)


def set_chem_addition_op_conditions(blk, **kwargs):
    m = blk.model()
    m.db.get_unit_operation_parameters("chemical_addition")
    blk.unit.load_parameters_from_database()


def add_costing(m):
    m.fs.costing = ZeroOrderCosting(
        case_study_definition="src/wrd/meta_data/wrd_case_study.yaml"
    )


def add_chem_addition_costing(m, blk, flowsheet_costing_block=None):
    if flowsheet_costing_block is None:
        flowsheet_costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block
    )


def init_chem_addition(blk, solver=None):
    if solver is None:
        solver = get_solver()

    # print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    # print(f"Chem Addition Degrees of Freedom: {degrees_of_freedom(m.fs.chem_addition)}")

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()

    propagate_state(blk.unit_to_product)
    blk.product.initialize()


def print_chem_addition_costing_breakdown(blk):
    print(
        f'{"Hydrogen Peroxide Dose":<35s}{f"{blk.unit.chemical_dosage[0]():<25,.0f} mg/L"}'
    )
    print(
        f'{"Chem Addition Capital Cost":<35s}{f"${blk.unit.costing.capital_cost():<25,.0f}"}'
    )


def solve(m, solver=None, tee=True, raise_on_failure=True):

    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)

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


def main():
    m = build_system()
    set_system_conditions(m.fs.chem_addition)
    set_chem_addition_op_conditions(m.fs.chem_addition)
    set_chem_addition_scaling(m.fs.chem_addition, calc_blk_scaling_factors=True)
    init_chem_addition(m.fs.chem_addition)
    solve(m)

    add_costing(m)
    add_chem_addition_costing(m, m.fs.chem_addition)

    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    solve(m)

    m.fs.costing.display()


if __name__ == "__main__":
    main()