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


# test
def build_decarbonator(blk, prop_package):

    blk.feed = Feed(property_package=prop_package)
    blk.product = Product(property_package=prop_package)

    blk.unit = StateJunction(property_package=prop_package)

    blk.unit.power_consumption = Var(
        initialize=0,
        domain=NonNegativeReals,
        units=pyunits.kW,
        doc="Power consumption of decarbonator",
    )
