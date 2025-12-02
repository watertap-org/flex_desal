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
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models import Product, Feed, StateJunction, Mixer, Separator
from idaes.core.util.model_statistics import *
from watertap.core.util.initialization import *
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.unit_models.zero_order.ultra_filtration_zo import UltraFiltrationZO
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.unit_models.pressure_changer import Pump

from watertap.core import Database


# TODO:
# 1. Unfix the variable energy_electric_flow_vol_inlet


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.db = Database(dbpath="watertap/flowsheets/flex_desal/wrd/meta_data")
    m.fs.properties = WaterParameterBlock(solute_list=["tds", "tss"])

    m.fs.UF = FlowsheetBlock(dynamic=False)
    build_UF(m.fs.UF, m.fs.properties)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def build_UF(blk, prop_package) -> None:

    print(f'\n{"=======> BUILDING ULTRAFILTRATION SYSTEM <=======":^60}\n')

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.UF_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)

    blk.unit = UltraFiltrationZO(property_package=prop_package, database=m.db)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_disposal = Arc(
        source=blk.unit.byproduct,
        destination=blk.disposal.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.treated,
        destination=blk.product.inlet,
    )


def set_inlet_conditions(blk):
    blk.feed.properties[0.0].flow_mass_comp["H2O"].fix(171.37)
    blk.feed.properties[0.0].flow_mass_comp["tds"].fix(600)
    blk.feed.properties[0.0].flow_mass_comp["tss"].fix(5.22e-6)


def set_UF_op_conditions(blk):
    # print(f"UF Degrees of Freedom: {degrees_of_freedom(blk)}")
    # blk.unit.recovery_frac_mass_H2O.fix(0.99)
    # blk.unit.removal_frac_mass_comp[0, "tds"].fix(1e-3)
    # blk.unit.removal_frac_mass_comp[0, "tss"].fix(0.9)
    # blk.unit.energy_electric_flow_vol_inlet.fix(0.05)
    load_parameters(blk)


def init_UF(blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    print(
        "\n\n-------------------- INITIALIZING ULTRAFILTRATION --------------------\n\n"
    )

    print(f"UF Degrees of Freedom: {degrees_of_freedom(blk)}")
    print("\n\n")

    # assert_no_degrees_of_freedom(m)

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize(optarg=optarg)
    propagate_state(blk.unit_to_disposal)
    propagate_state(blk.unit_to_product)

    blk.product.initialize(optarg=optarg)
    blk.disposal.initialize(optarg=optarg)


def add_UF_costing(m, blk, costing_blk=None):
    if costing_blk is None:
        costing_blk = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_blk)


def add_UF_scaling(blk):
    set_scaling_factor(blk.disposal.properties[0.0].flow_mass_comp["tds"], 1e3)
    set_scaling_factor(blk.unit.properties_byproduct[0.0].flow_mass_comp["tds"], 1e3)


def load_parameters(blk):
    m = blk.model()
    m.db.get_unit_operation_parameters("ultra_filtration")
    blk.unit.load_parameters_from_database()


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
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = build_system()
    set_inlet_conditions(m.fs.UF)
    set_UF_op_conditions(m.fs.UF)
    add_UF_scaling(m.fs.UF)

    init_UF(m.fs.UF)
    solve(m)

    add_costing(m)
    add_UF_costing(m, m.fs.UF)
    m.fs.costing.cost_process()
    m.fs.costing.initialize()

    solve(m)

    print("Degrees of Freedom: ", degrees_of_freedom(m))
