import os
import yaml
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    Param,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.models.unit_models import Product, Feed, StateJunction
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.solvers import get_solver
from watertap.costing import WaterTAPCosting
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock

from models import ChemicalAddition
from srp.utils import touch_flow_and_conc
from wrd.utilities import get_config_value, load_config, get_config_file

__all__ = [
    "build_chem_addition",
    "set_chem_addition_op_conditions",
    "add_chem_addition_costing",
    "report_chem_addition",
    "initialize_chem_addition",
]

solver = get_solver()


def get_chem_data(chem_data, chemical_name, default=None):
    chem_config = {}

    chem_config["ratio_in_solution"] = get_config_value(
        chem_data,
        "ratio_in_solution",
        chemical_name,
    )

    chem_config["solution_density"] = get_config_value(
        chem_data,
        "solution_density",
        chemical_name,
    )

    chem_config["unit_cost"] = get_config_value(
        chem_data,
        "unit_cost",
        chemical_name,
    )

    chem_config["chemical_dosage"] = get_config_value(
        chem_data,
        "chemical_dosage",
        chemical_name,
    )

    return chem_config


def build_system(chemical_name=None):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2021
    m.fs.costing.base_period = pyunits.month
    m.fs.properties = NaClParameterBlock()

    config_file_name = "chemical_addition.yaml"
    config = get_config_file(config_file_name)
    m.fs.chem_data = load_config(config)

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


def build_chem_addition(blk, chemical_name=None, prop_package=None, file=None):

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.properties

    if chemical_name is None:
        chemical_name = "default_chemical"

    name = chemical_name.replace("_", " ").upper()
    print(f'\n{f"=======> BUILDING {name} ADDITION UNIT <=======":^60}\n')

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.chem_config = get_chem_data(m.fs.chem_data, chemical_name, None)

    blk.unit = ChemicalAddition(
        property_package=prop_package,
        chemical=chemical_name,
        chemical_data=blk.chem_config,
    )

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


def set_inlet_conditions(m, Qin=2637, Cin=0.5):

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
        dose = blk.chem_config["chemical_dosage"]
        if dose is None:
            raise ValueError("dose must be provided to set_chem_addition_op_conditions")

    blk.unit.dose.fix(dose)


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


def add_chem_addition_costing(
    blk, costing_package=None, chem_cost=None, chem_purity=None
):

    if chem_cost is None:
        chem_cost = blk.chem_config["unit_cost"]
        if chem_cost is None:
            raise ValueError("chem_cost must be provided to add_chem_addition_costing")

    if chem_purity is None:
        chem_purity = blk.chem_config["ratio_in_solution"]
        if chem_purity is None:
            chem_purity = 1.0  # assume 100% purity if not provided

    m = blk.model()
    if costing_package is None:
        costing_package = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_package)

    cost_var = getattr(m.fs.costing, f"{blk.unit.config.chemical}")
    cost_var.cost.fix(chem_cost)


def report_chem_addition(blk, w=35):
    chem_name = blk.unit.config.chemical.replace("_", " ").title()
    feed_flow = pyunits.convert(
        blk.unit.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.min,
    )
    title = f"{chem_name} Addition Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(f'{f"Inlet Flow":<{w}s}{f"{value(feed_flow):<{w}.2f}gpm"}')
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
        f'{f"{chem_name} Pump":<{w}s}{value(blk.unit.pumping_power):<{w}.3e}{f"{pyunits.get_units(blk.unit.pumping_power)}"}'
    )

    m = blk.model()
    cost_var = getattr(m.fs.costing, f"{blk.unit.config.chemical}")
    unit_cost = pyunits.convert(cost_var.cost, to_units=pyunits.USD_2021 / pyunits.kg)
    print(
        f'{f"{chem_name} unit cost":<{w}s}{value(unit_cost):<{w}.3f}{f"{pyunits.get_units(unit_cost)}"}'
    )
    print(
        f'{"Chem Addition Operating Cost":<{w}s}{f"${value(m.fs.costing.aggregate_flow_costs[blk.unit.config.chemical]):<{w}.3f}$/month"}'
    )


def main(
    chemical_name="ammonium_sulfate",
    Qin=2637,
    Cin=0.5,
    # If hard coding, need to pass units somewhere
    dose=None,  # 0.01,
    chem_cost=None,  # 0.5,
    chem_purity=None,  # 0.9,
):
    m = build_system(chemical_name=chemical_name)
    # Add units to chem_cost after costing system defines currency units
    if chem_cost is not None:
        chem_cost = chem_cost * pyunits.USD_2021 / pyunits.kg
    add_chem_addition_costing(
        m.fs.chem_addition, chem_cost=chem_cost, chem_purity=chem_purity
    )
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin)
    set_chem_addition_op_conditions(m.fs.chem_addition, dose=dose)
    initialize_system(m)
    m.fs.costing.cost_process()
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_chem_addition(m.fs.chem_addition, w=40)
    return m


if __name__ == "__main__":
    chem = "scale_inhibitor"
    m = main(chemical_name=chem)
