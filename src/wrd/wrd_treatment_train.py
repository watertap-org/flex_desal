from os import path
from pyomo.environ import (
    ConcreteModel,
    Param,
    check_optimal_termination,
    value,
    assert_optimal_termination,
    units as pyunits,
    value,
    TransformationFactory,
)

import pyomo.contrib.parmest.parmest as parmest
from pyomo.contrib.parmest.experiment import Experiment

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    MixingType,
    MomentumMixingType,
    Mixer,
    Separator,
    Product,
    Feed,
)

from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock

from wrd.components.UF import *
from wrd.components.chemical_addition import *
from wrd.components.translator_ZO_to_NaCl import (
    TranslatorZOtoNaCl,
)
from wrd.components.ro_system import *
from wrd.components.decarbonator import *
from wrd.components.UV_aop import *


def build_wrd_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Get working directory path
    dir_path = os.path.dirname(os.path.abspath(__file__))
    m.db = Database(dbpath=os.path.join(dir_path, "meta_data"))

    m.fs.properties = WaterParameterBlock(solute_list=["tds", "tss"])

    # RO properties
    m.fs.ro_properties = NaClParameterBlock()

    # Add units
    m.fs.feed = Feed(property_package=m.fs.properties)

    # Chemical addition units
    m.fs.ammonia_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(m.fs.ammonia_addition, "ammonia", m.fs.properties)

    m.fs.hypochlorite_addition = FlowsheetBlock(dynamic=False)
    build_chem_addition(
        m.fs.hypochlorite_addition, "sodium_hypochlorite", m.fs.properties
    )

    # UF unit
    m.fs.UF = FlowsheetBlock(dynamic=False)
    build_UF(m.fs.UF, m.fs.properties)

    # Translator block between ZO to RO property packages
    m.fs.translator_ZO_to_RO = TranslatorZOtoNaCl(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.ro_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # RO unit
    m.fs.ro_system = FlowsheetBlock(dynamic=False)
    build_wrd_ro_system(m.fs.ro_system, prop_package=m.fs.ro_properties)

    # UV AOP - Still using ro_properties
    m.fs.UV_aop = FlowsheetBlock(dynamic=False)
    build_UV(m.fs.UV_aop, prop_package=m.fs.ro_properties)

    # Decarbonator
    m.fs.decarbonator = FlowsheetBlock(dynamic=False)
    build_decarbonator(m.fs.decarbonator)

    m.fs.product = Product(property_package=m.fs.properties)

    return m


def add_connections(m):
    # Connect feed to ammonia addition
    m.fs.feed_to_chem_add = Arc(
        source=m.fs.feed.outlet, destination=m.fs.ammonia_addition.feed.inlet
    )
    # Connect ammonia addition to hypochlorite addition
    m.fs.ammonia_to_hypochlorite = Arc(
        source=m.fs.ammonia_addition.product.outlet,
        destination=m.fs.hypochlorite_addition.feed.inlet,
    )
    # Connect hypochlorite addition to UF
    m.fs.hypochlorite_to_UF = Arc(
        source=m.fs.hypochlorite_addition.product.outlet, destination=m.fs.UF.feed.inlet
    )
    # Connect UF to RO translator
    m.fs.UF_to_translator = Arc(
        source=m.fs.UF.product.outlet, destination=m.fs.translator_ZO_to_RO.inlet
    )
    # Connect RO translator to RO
    m.fs.translator_to_ro = Arc(
        source=m.fs.translator_ZO_to_RO.outlet, destination=m.fs.ro_system.feed.inlet
    )
    # Connect RO to UV_aop
    m.fs.ro_to_uv = Arc(
        source=m.fs.ro_system.product.outlet, destination=m.fs.UV_aop.feed.inlet
    )
    # Connect UV_aop to Decarbonator
    m.fs.uv_to_decarbonator = Arc(
        source=m.fs.UV_aop.product.outlet, destination=m.fs.decarbonator.feed.inlet
    )
    # Connect Decarbonator to the Product
    m.fs.decarbonator_to_product = Arc(
        source=m.fs.decarbonator.product.outlet, destination=m.fs.product.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_wrd_inlet_conditions(m):
    # Inlet conditions
    # m.fs.feed.properties[0].pressure.fix(101325)  # Fix feed pressure to 1 atm
    # m.fs.feed.properties[0].temperature.fix(298.15)  # Fix feed temperature to 25 C
    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(174)  # Fix feed water flow rate
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(0.05)  # Fix feed salt flow rate
    m.fs.feed.properties[0].flow_mass_comp["tss"].fix(0.1)  # Fix feed salt flow rate


def set_wrd_operating_conditions(m):
    # Operating conditions
    set_chem_addition_op_conditions(blk=m.fs.ammonia_addition)
    set_chem_addition_op_conditions(blk=m.fs.hypochlorite_addition)
    set_UF_op_conditions(m.fs.UF)
    set_ro_system_op_conditions(m.fs.ro_system)
    set_UV_aop_op_conditions(m.fs.UV_aop)
    set_decarbonator_op_conditions(m.fs.decarbonator)


def initialize_wrd_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_chem_add)
    init_chem_addition(m.fs.ammonia_addition)
    propagate_state(m.fs.ammonia_to_hypochlorite)
    init_chem_addition(m.fs.hypochlorite_addition)
    propagate_state(m.fs.hypochlorite_to_UF)
    init_UF(m.fs.UF)
    propagate_state(m.fs.s0UF_to_translator4)
    m.fs.translator_ZO_to_RO.initialize()
    propagate_state(m.fs.translator_to_ro)
    initialize_ro_system(m.fs.ro_system)
    propagate_state(m.fs.ro_to_uv)
    initialize_UV_aop(m.fs.UV_aop)
    propagate_state(m.fs.uv_to_decarbonator)
    initialize_decarbonator(m.fs.decarbonator)
    propagate_state(m.fs.decarbonator_to_product)
    m.fs.product.initialize()


def set_wrd_system_scaling(m):
    # Properties Scaling
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    # Does ZO property block also require scaling?

    set_chem_addition_scaling(blk=m.fs.ammonia_addition)
    set_chem_addition_scaling(blk=m.fs.hypochlorite_addition)
    add_UF_scaling(m.fs.UF)
    add_ro_scaling(m.fs.ro_system)
    add_UV_scaling(m.fs.UV_aop)
    add_decarbonator_scaling(m.fs.decarbonator)


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


if __name__ == "__main__":

    m = build_wrd_system()
    add_connections(m)
    set_wrd_inlet_conditions(m)
    set_wrd_operating_conditions(m)
    set_wrd_system_scaling(m)
    calculate_scaling_factors(m)
    initialize_wrd_system(m)
    m.fs.ro_system.feed.display()

    # assert False

    try:
        results = solve(m)
        assert_optimal_termination(results)
    except:
        print_infeasible_constraints(m)
