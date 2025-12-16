import os
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
from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock

from models import ChemicalAddition
from wrd.utilities import load_config, get_config_value, get_config_file
from srp.utils import touch_flow_and_conc

__all__ = [
    "build_chem_addition",
    "set_chem_addition_op_conditions",
    "add_chem_addition_costing",
    "report_chem_addition",
    "initialize_chem_addition",
]

solver = get_solver()


def build_system(chemical_name=None):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2007
    current_script_path = os.path.abspath(__file__)
    current_directory = os.path.dirname(current_script_path)
    parent_directory = os.path.dirname(current_directory)
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.chem_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m.fs.chem_addition, chemical_name, m.fs.properties)

    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.feed_to_chem_addition = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.chem_addition.feed.inlet,
    )

    m.fs.chem_addition_to_product = Arc(
        source=m.fs.chem_addition.product.outlet,
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


def build_chem_addition(blk, chemical_name=None, prop_package=None):

    if chemical_name is None:
        # raise ValueError("chemical_name must be provided to build_chem_addition")
        name = "default_chemical"
    else:
        name = chemical_name.replace("_", " ").upper()

    print(f'\n{f"=======> BUILDING {name} ADDITION UNIT <=======":^60}\n')

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = ChemicalAddition(property_package=prop_package, chemical=chemical_name)
    blk.product = StateJunction(property_package=prop_package)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )

    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_inlet_conditions(m, Qin=2637, Cin=0.5, file="wrd_ro_inputs_8_19_21.yaml"):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + 27,
        },
        hold_state=True,
    )


def set_chem_addition_op_conditions(blk, dose=None):

    if dose is None:
        raise ValueError("dose must be provided to set_chem_addition_op_conditions")

    blk.unit.dose.fix(dose)


def add_chem_addition_costing(blk, flowsheet_costing_block=None, cost_capital=False):

    if flowsheet_costing_block is None:
        m = blk.model()
        flowsheet_costing_block = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block,
        costing_method_arguments={"cost_capital": cost_capital},
    )


def report_chem_addition(blk, w=30):
    chem_name = blk.unit.config.chemical.replace("_", " ").title()
    title = f"{chem_name} Addition Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{f"{chem_name} Dose":<{w}s}{f"{value(pyunits.convert(blk.unit.dose, to_units=pyunits.mg/pyunits.liter)):<{w}.2f}mg/L"}'
    )
    print(
        f'{f"{chem_name} Mass Flow":<{w}s}{value(blk.unit.chemical_flow_mass):<{w}.3e}{f"{pyunits.get_units(blk.unit.chemical_flow_mass)}"}'
    )
    print(
        f'{f"{chem_name} Vol. Flow":<{w}s}{value(blk.unit.chemical_soln_flow_vol):<{w}.3e}{f"{pyunits.get_units(blk.unit.chemical_soln_flow_vol)}"}'
    )
    print(
        f'{f"{chem_name} Pumping Power":<{w}s}{value(blk.unit.pumping_power):<{w}.3e}{f"{pyunits.get_units(blk.unit.pumping_power)}"}'
    )
    # print(
    #     f'{"Chem Addition Capital Cost":<{w}s}{f"${blk.unit.costing.capital_cost():<{w}.3f}"}'
    # )


def initialize_system(blk):
    blk.fs.feed.initialize()
    propagate_state(blk.fs.feed_to_chem_addition)

    initialize_chem_addition(blk.fs.chem_addition)

    propagate_state(blk.fs.chem_addition_to_product)
    blk.fs.product.initialize()


def initialize_chem_addition(blk):
    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    blk.unit.initialize()

    propagate_state(blk.unit_to_product)
    blk.product.initialize()


def set_system_operating_conditions(m, Qin=2637, Cin=0.5, dose=0.01):
    set_inlet_conditions(m, Qin=Qin, Cin=Cin)
    set_chem_addition_op_conditions(m.fs.chem_addition, dose=dose)


def main():

    m = build_system()
    add_chem_addition_costing(m.fs.chem_addition)
    calculate_scaling_factors(m)
    set_system_operating_conditions(m, Qin=2637, Cin=0.5, dose=0.010)
    initialize_system(m)
    m.fs.costing.cost_process()

    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_chem_addition(m.fs.chem_addition, w=40)

    return m


if __name__ == "__main__":
    m = main()
