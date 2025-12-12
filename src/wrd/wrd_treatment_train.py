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
from wrd.components.translator_ZO_to_NaCl import TranslatorZOtoNaCl
from wrd.components.translator_NaCl_to_ZO import TranslatorNaCltoZO
from wrd.components.ro_system import *
from wrd.components.decarbonator import *
from wrd.components.uv_aop import *
from wrd.utilities import load_config, get_config_file, get_config_value
from srp.components.generic_separator import *


def build_wrd_system(number_stages=3, **kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Get working directory path
    dir_path = os.path.dirname(os.path.abspath(__file__))
    m.db = Database(dbpath=os.path.join(dir_path, "meta_data"))

    config_file_name = get_config_file("wrd_feed_flow.yaml")
    m.fs.config_data = load_config(config_file_name)

    # ZO Properties
    m.fs.properties = WaterParameterBlock(solute_list=["tds"])
    # RO properties
    m.fs.ro_properties = NaClParameterBlock()

    # Add units
    m.fs.feed = Feed(property_package=m.fs.properties)

    # Pre- UF Treatment chemical addition units (read from metadata)
    m.fs.pre_treat_chem_list = [
        "ammonium_sulfate",
        "sodium_hypochlorite",
        "sulfuric_acid",
        "scale_inhibitor",
    ]
    for chem_name in m.fs.pre_treat_chem_list:
        m.fs.add_component(chem_name + "_addition", FlowsheetBlock(dynamic=False))
        build_chem_addition(
            m.fs.find_component(chem_name + "_addition"), chem_name, m.fs.properties
        )

    # UF unit
    m.fs.UF = FlowsheetBlock(dynamic=False)
    # build_UF(m.fs.UF, m.fs.properties)
    build_separator(
        blk=m.fs.UF, prop_package=m.fs.ro_properties, outlet_list=["to_RO", "to_waste"]
    )

    # Translator block between ZO to RO property packages
    m.fs.translator_ZO_to_RO = TranslatorZOtoNaCl(
        inlet_property_package=m.fs.properties,
        outlet_property_package=m.fs.ro_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # RO unit
    m.fs.ro_system = FlowsheetBlock(dynamic=False)
    number_stages = 3
    if "number_stages" in kwargs:
        number_stages = kwargs["number_stages"]
    build_wrd_ro_system(
        m.fs.ro_system,
        prop_package=m.fs.ro_properties,
        number_stages=number_stages,
    )

    # UV AOP - Still using ro_properties
    m.fs.UV_aop = FlowsheetBlock(dynamic=False)
    build_uv_aop(m.fs.UV_aop, prop_package=m.fs.ro_properties)

    # Decarbonator - Still using ro_properties
    m.fs.decarbonator = FlowsheetBlock(dynamic=False)
    build_decarbonator(m.fs.decarbonator, prop_package=m.fs.ro_properties)

    m.fs.translator_RO_to_ZO = TranslatorNaCltoZO(
        inlet_property_package=m.fs.ro_properties,
        outlet_property_package=m.fs.properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # Post-Treatment chemical addition units - ZO Models
    m.fs.post_treat_chem_list = [
        "calcium_hydroxide",
        "sodium_hydroxide",
        "sodium_hypochlorite_post",
    ]

    for chem_name in m.fs.post_treat_chem_list:
        m.fs.add_component(chem_name + "_addition", FlowsheetBlock(dynamic=False))
        build_chem_addition(
            m.fs.find_component(chem_name + "_addition"), chem_name, m.fs.properties
        )

    m.fs.product = Product(property_package=m.fs.properties)

    # Combined chemical list for operating conditions, scaling, and costing(?)
    m.fs.chemical_list = list(m.fs.pre_treat_chem_list) + list(
        m.fs.post_treat_chem_list
    )
    return m


def add_wrd_connections(m):
    # Connect pre-UF chemical chain: feed -> chem1 -> chem2 -> ... -> UF
    for i in range(len(m.fs.pre_treat_chem_list)):
        chem_name = m.fs.pre_treat_chem_list[i]
        unit = m.fs.find_component(chem_name + "_addition")
        if i == 0:
            # Connect feed to first chemical
            m.fs.add_component(
                f"feed_to_{chem_name}",
                Arc(
                    source=m.fs.feed.outlet,
                    destination=unit.feed.inlet,
                ),
            )
        else:
            # Connect each chemical to the next
            prev = m.fs.pre_treat_chem_list[i - 1]
            prev_unit = m.fs.find_component(prev + "_addition")
            m.fs.add_component(
                f"{prev}_to_{chem_name}",
                Arc(
                    source=prev_unit.product.outlet,
                    destination=unit.feed.inlet,
                ),
            )

    # Connect last pre treat chemical to UF_pump
    m.fs.add_component(
        f"{m.fs.pre_treat_chem_list[-1]}_to_translator",
        Arc(
            source=m.fs.find_component(
                m.fs.pre_treat_chem_list[-1] + "_addition"
            ).product.outlet,
            destination=m.fs.translator_ZO_to_RO.inlet,
        ),
    )

    # Connect RO translator to RO
    m.fs.translator_to_uf = Arc(
        source=m.fs.translator_ZO_to_RO.outlet, destination=m.fs.UF.feed.inlet
    )
    # Connect UF Pump to UF

    # UF to RO translator
    m.fs.UF_to_ro = Arc(
        source=m.fs.UF.to_RO.outlet, destination=m.fs.ro_system.feed.inlet
    )

    # Connect RO to UV_aop
    m.fs.ro_to_uv = Arc(
        source=m.fs.ro_system.permeate.outlet, destination=m.fs.UV_aop.feed.inlet
    )

    # Connect UV_aop to Decarbonator
    m.fs.uv_to_decarbonator = Arc(
        source=m.fs.UV_aop.product.outlet, destination=m.fs.decarbonator.feed.inlet
    )

    # Connect Decarbonator to translator
    m.fs.ro_to_translator = Arc(
        source=m.fs.decarbonator.product.outlet,
        destination=m.fs.translator_RO_to_ZO.inlet,
    )

    # Chain post-treatment chemicals (decarbonator -> post1 -> post2 -> ... -> product)
    for i in range(len(m.fs.post_treat_chem_list)):
        chem_name = m.fs.post_treat_chem_list[i]
        if i == 0:
            # Connect decarb to first chemical
            m.fs.add_component(
                "decarb_to_" + chem_name,
                Arc(
                    source=m.fs.translator_RO_to_ZO.outlet,
                    destination=m.fs.find_component(chem_name + "_addition").feed.inlet,
                ),
            )
        else:  #
            # Connect each chemical to the next
            m.fs.add_component(
                m.fs.post_treat_chem_list[i - 1] + "_to_" + chem_name,
                Arc(
                    source=m.fs.find_component(
                        m.fs.post_treat_chem_list[i - 1] + "_addition"
                    ).product.outlet,
                    destination=m.fs.find_component(chem_name + "_addition").feed.inlet,
                ),
            )

    # Connect last chemical to the Product
    m.fs.add_component(
        m.fs.post_treat_chem_list[-1] + "_to_product",
        Arc(
            source=m.fs.find_component(
                m.fs.post_treat_chem_list[-1] + "_addition"
            ).product.outlet,
            destination=m.fs.product.inlet,
        ),
    )
    TransformationFactory("network.expand_arcs").apply_to(m)


def set_wrd_inlet_conditions(m):
    # Inlet conditions
    number_trains = m.fs.ro_system.number_trains
    Qin = get_config_value(
        m.fs.config_data,
        "feed_flow_water",
        "feed_stream",
    )

    Cin = get_config_value(
        m.fs.config_data, "feed_conductivity", "feed_stream"
    ) * get_config_value(
        m.fs.config_data, "feed_conductivity_conversion", "feed_stream"
    )

    # Would like to load Pin and Tin as well once property model is changed.

    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(
        number_trains * feed_mass_flow_water
    )  # Fix feed water flow rate
    # Not sure this is correct way to translate salinate to tds and tss
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(
        feed_mass_flow_salt
    )  # Fix feed salt flow rate
    # m.fs.feed.properties[0].flow_mass_comp["tss"].fix(0)  # Fix feed salt flow rate


def set_wrd_operating_conditions(m):
    # Operating conditions
    for chem_name in m.fs.chemical_list:
        set_chem_addition_op_conditions(
            blk=m.fs.find_component(chem_name + "_addition")
        )
    # set_UF_op_conditions(m.fs.UF)
    uf_splits = {
        "to_RO": {"H2O": 0.99, "NaCl": 0.99},
    }
    set_separator_op_conditions(m.fs.UF, split_fractions=uf_splits)
    set_ro_system_op_conditions(m.fs.ro_system)
    set_uv_aop_op_conditions(m.fs.UV_aop)
    set_decarbonator_op_conditions(m.fs.decarbonator)


def initialize_wrd_system(m):
    m.fs.feed.initialize()
    # Initialize pre-UF chemical chain
    for i, chem_name in enumerate(m.fs.pre_treat_chem_list):
        if i == 0:
            propagate_state(m.fs.find_component("feed_to_" + chem_name))
        else:
            prev = m.fs.pre_treat_chem_list[i - 1]
            propagate_state(m.fs.find_component(prev + "_to_" + chem_name))

        init_chem_addition(m.fs.find_component(chem_name + "_addition"))

    # propagate from last pre-UF chemical to UF
    propagate_state(
        m.fs.find_component(m.fs.pre_treat_chem_list[-1] + "_to_translator")
    )
    # init_UF(m.fs.UF)
    m.fs.translator_ZO_to_RO.initialize()
    propagate_state(m.fs.translator_to_uf)
    init_separator(m.fs.UF)
    propagate_state(m.fs.UF_to_ro)

    # propagate last pre-RO to translator
    # propagate_state(
    #     m.fs.find_component(m.fs.pre_treat_chem_list[-1] + "_to_translator")
    # )

    # m.fs.translator_ZO_to_RO.initialize()
    # propagate_state(m.fs.translator_to_ro)
    initialize_ro_system(m.fs.ro_system)
    propagate_state(m.fs.ro_to_uv)
    initialize_uv_aop(m.fs.UV_aop)
    propagate_state(m.fs.uv_to_decarbonator)
    initialize_decarbonator(m.fs.decarbonator)

    # Initialize post-treatment chemical chain (downstream of decarbonator)
    for i, chem_name in enumerate(m.fs.post_treat_chem_list):
        if i == 0:
            propagate_state(m.fs.find_component("decarb_to_" + chem_name))
        else:
            prev = m.fs.post_treat_chem_list[i - 1]
            propagate_state(m.fs.find_component(prev + "_to_" + chem_name))
        init_chem_addition(m.fs.find_component(chem_name + "_addition"))

    propagate_state(m.fs.find_component(m.fs.post_treat_chem_list[-1] + "_to_product"))
    m.fs.product.initialize()


def set_wrd_system_scaling(m):
    # Properties Scaling
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    # Does ZO property block also require scaling?
    for chem_name in m.fs.chemical_list:
        set_chem_addition_scaling(blk=m.fs.find_component(chem_name + "_addition"))

    # add_UF_scaling(m.fs.UF)
    add_ro_scaling(m.fs.ro_system)
    add_uv_aop_scaling(m.fs.UV_aop)
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


def main(number_stages=3, date="8_19_21"):
    m = build_wrd_system(number_stages=number_stages, date=date)
    add_wrd_connections(m)
    set_wrd_inlet_conditions(m)
    set_wrd_operating_conditions(m)
    set_wrd_system_scaling(m)
    calculate_scaling_factors(m)
    initialize_wrd_system(m)
    try:
        results = solve(m)
        assert_optimal_termination(results)
    except:
        print_infeasible_constraints(m)
    return m


if __name__ == "__main__":
    number_stages = 3
    date = "8_19_21"
    main(number_stages=number_stages, date=date)
    # m = build_wrd_system(number_stages=number_stages, date=date)
    # add_connections(m)
    # set_wrd_inlet_conditions(m)
    # set_wrd_operating_conditions(m)
    # set_wrd_system_scaling(m)
    # calculate_scaling_factors(m)
    # initialize_wrd_system(m)
    # m.fs.ro_system.feed.display()

    # # assert False

    # try:
    #     results = solve(m)
    #     assert_optimal_termination(results)
    # except:
    #     print_infeasible_constraints(m)

    # # report_wrd(m)
