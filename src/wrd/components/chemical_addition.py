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

__all__ = [
    "build_chem_addition",
    "set_chem_addition_op_conditions",
    "add_chem_addition_costing",
    "report_chem_addition",
    "initialize_chem_addition",
]

solver = get_solver()


current_script_path = os.path.abspath(__file__)
current_directory = os.path.dirname(current_script_path)
parent_directory = os.path.dirname(current_directory)
default_chem_addition_config_file = os.path.join(
    parent_directory, "meta_data", "chemical_addition.yaml"
)


def build_system(chemical_name=None):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
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


def build_chem_addition(blk, chemical_name=None, prop_package=None, file=None):

    if chemical_name is None:
        # raise ValueError("chemical_name must be provided to build_chem_addition")
        name = "default_chemical"
    else:
        name = chemical_name.replace("_", " ").upper()

    if file is None:
        file = default_chem_addition_config_file

    with open(file, "r") as f:
        data = yaml.safe_load(f)

    blk.chem_data = data.get(chemical_name, {})

    config_data = {}

    if blk.chem_data != {}:
        config_data["solution_density"] = blk.chem_data["solution_density"]["value"]
        config_data["ratio_in_solution"] = blk.chem_data["ratio_in_solution"]["value"]

    print(f'\n{f"=======> BUILDING {name} ADDITION UNIT <=======":^60}\n')

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.properties

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = ChemicalAddition(
        property_package=prop_package,
        chemical=chemical_name,
        chemical_data=config_data,
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
        dose = blk.chem_data.get("chemical_dosage", None)
        if dose is not None:
            dose = dose.get("value", None)
        elif dose is None:
            raise ValueError(
                "dose must be provided to set_chem_addition_op_conditions or in the config file"
            )

    blk.unit.dose.fix(dose)


def report_chem_addition(blk, w=30):
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


def add_chem_addition_costing(
    blk, costing_package=None, chem_cost=None, chem_purity=None, chem_cost_units=None
):
    if chem_cost is None:
        chem_cost = blk.chem_data["unit_cost"]["value"]
        if chem_cost is None:
            raise ValueError("chem_cost must be provided to add_chem_addition_costing")

    if chem_purity is None:
        chem_purity = blk.chem_data.get("purity", None)
        if chem_purity is None:
            chem_purity = 1.0  # assume 100% purity if not provided

    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    if chem_cost_units is None:
        try:
            chem_cost_units = blk.chem_data["unit_cost"]["units"]
            num = getattr(pyunits, chem_cost_units.split("/")[0])
            den = getattr(pyunits, chem_cost_units.split("/")[1])
            chem_cost_units = num / den
        except:
            chem_cost_units = costing_package.base_currency / pyunits.kg

    if blk.unit.config.chemical not in costing_package._registered_flows.keys():
        blk.unit.cost = Param(
            initialize=chem_cost,
            units=chem_cost_units,
            mutable=True,
            doc=f"{blk.unit.config.chemical.replace('_', ' ').title()} cost",
        )
        blk.unit.purity = Param(
            initialize=chem_purity,
            units=pyunits.dimensionless,
            mutable=True,
            doc=f"{blk.unit.config.chemical.replace('_', ' ').title()} purity",
        )
        costing_package.register_flow_type(
            blk.unit.config.chemical, blk.unit.cost / blk.unit.purity
        )

    costing_package.cost_flow(blk.unit.chemical_flow_mass, blk.unit.config.chemical)
    costing_package.cost_flow(blk.unit.pumping_power, "electricity")


def main(
    chemical_name="ammonia",
    Qin=2637,
    Cin=0.5,
    dose=0.01,
    chem_cost=0.5,
    chem_purity=0.9,
):

    m = build_system(chemical_name=chemical_name)
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
    chem = "caustic"
    m = main(chemical_name=chem)

    # import yaml
    # _f = "/Users/ksitterl/Documents/Python/flex_desal/flex_desal/src/wrd/meta_data/chemical_addition.yaml"
    # with open(_f, "r") as f:
    #     data = yaml.safe_load(f)

    # print(data)
    # print(type(data))  # <class 'dict'>
    # print(data[m.fs.chem_addition.unit.config.chemical])
    # print(parent_directory)
    # print(default_chem_addition_config_file)
