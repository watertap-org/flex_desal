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

# from watertap.flowsheets.flex_desal.wrd.components.ro_system import *
from wrd.components.ro_system import *
from wrd.components.translator_ZO_to_NaCl import (
    TranslatorZOtoNaCl,
)


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
    m.fs.product = Product(property_package=m.fs.properties)

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
    m.fs.ro_train = FlowsheetBlock(dynamic=False)
    build_wrd_ro_system(
        m.fs.ro_train, prop_package=m.fs.ro_properties
    )  # Need to pass stage_num now

    # Translator block between RO to ZO property packages

    return m


def add_connections(m):
    # Connect feed to ammonia addition
    m.fs.s01 = Arc(
        source=m.fs.feed.outlet, destination=m.fs.ammonia_addition.feed.inlet
    )
    # Connect ammonia addition to hypochlorite addition
    m.fs.s02 = Arc(
        source=m.fs.ammonia_addition.product.outlet,
        destination=m.fs.hypochlorite_addition.feed.inlet,
    )
    # Connect hypochlorite addition to UF
    m.fs.s03 = Arc(
        source=m.fs.hypochlorite_addition.product.outlet, destination=m.fs.UF.feed.inlet
    )
    # Connect UF to RO translator
    m.fs.s04 = Arc(
        source=m.fs.UF.product.outlet, destination=m.fs.translator_ZO_to_RO.inlet
    )
    # Connect RO translator to RO
    m.fs.s05 = Arc(
        source=m.fs.translator_ZO_to_RO.outlet, destination=m.fs.ro_train.feed.inlet
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
    # set_ro_operation_conditions(m.fs.ro_train)
    set_ro_system_op_conditions(m.fs.ro_train)


def initialize_wrd_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.s01)
    init_chem_addition(m.fs.ammonia_addition)
    propagate_state(m.fs.s02)
    init_chem_addition(m.fs.hypochlorite_addition)
    propagate_state(m.fs.s03)
    init_UF(m.fs.UF)
    propagate_state(m.fs.s04)
    m.fs.translator_ZO_to_RO.initialize()
    propagate_state(m.fs.s05)
    initialize_ro_system(m.fs.ro_train)
    # m.fs.ro_train.total_ro_feed.initialize()
    # build_ro_inlet_stream(m.fs.ro_train, test=False)
    # initialize_ro_units(m.fs.ro_train)


def set_wrd_system_scaling(m):

    set_chem_addition_scaling(blk=m.fs.ammonia_addition)
    set_chem_addition_scaling(blk=m.fs.hypochlorite_addition)
    add_UF_scaling(m.fs.UF)
    add_ro_scaling(m.fs.ro_train)


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
    initialize_wrd_system(m)

    # m.fs.ro_train.total_ro_feed.display()
    m.fs.ro_train.feed.display()

    # assert False

    try:
        results = solve(m)
        assert_optimal_termination(results)
    except:
        print_infeasible_constraints(m)
