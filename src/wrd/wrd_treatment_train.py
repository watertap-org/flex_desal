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
from idaes.models.unit_models import Product, Feed
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.zero_order_properties import WaterParameterBlock

from wrd.components.chemical_addition import *
from wrd.components.ro_system import *

# from wrd.components.ro_system_new import build_ro_system
from wrd.components.decarbonator import *
from wrd.components.uv_aop import *

# from wrd.components.UF_feed_pumps import *
from wrd.components.pump import *
from wrd.components.UF_system import *
from wrd.components.ro_system import *
from wrd.components.ro_stage import *
from wrd.components.chemical_addition import *
from wrd.utilities import load_config, get_config_file, get_config_value
from srp.utils import touch_flow_and_conc


def build_wrd_system(num_pro_trains=4, num_tsro_trains=None, num_stages=2):

    if num_tsro_trains is None:
        num_tsro_trains = num_pro_trains

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    config_file_name = get_config_file("wrd_feed_flow.yaml")
    m.fs.config_data = load_config(config_file_name)
    m.fs.properties = NaClParameterBlock()

    # Add units
    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

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

    # UF
    build_uf_system(m=m, num_trains=num_pro_trains, prop_package=m.fs.properties)

    # PRO System
    build_ro_system(
        m=m,
        num_trains=num_pro_trains,
        num_stages=num_stages,
        prop_package=m.fs.properties,
    )

    # TSRO System
    m.fs.tsro_trains = Set(initialize=range(1, num_tsro_trains + 1))
    m.fs.tsro_header = StateJunction(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.tsro_header)
    m.fs.tsro_feed_separator = Separator(
        property_package=m.fs.properties,
        outlet_list=[f"to_tsro{i}" for i in m.fs.tsro_trains],
        split_basis=SplittingType.componentFlow,
    )
    m.fs.tsro_feed_separator.even_split = 1 / len(m.fs.tsro_trains)
    touch_flow_and_conc(m.fs.tsro_feed_separator)
    m.fs.tsro_train = FlowsheetBlock(m.fs.tsro_trains, dynamic=False)
    for t in m.fs.tsro_trains:
        build_ro_stage(m.fs.tsro_train[t], stage_num=3, prop_package=m.fs.properties)

    m.fs.tsro_product_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=[f"tsro{i}_to_product" for i in m.fs.tsro_trains],
        momentum_mixing_type=MomentumMixingType.none,
    )

    ro_system_prod_mixer_inlet_list = ["from_pro_product"] + [
        f"tsro{i}_to_ro_product" for i in m.fs.tsro_trains
    ]

    m.fs.ro_system_product_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=ro_system_prod_mixer_inlet_list,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.tsro_brine_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=[f"tsro{i}_to_brine" for i in m.fs.tsro_trains],
        momentum_mixing_type=MomentumMixingType.none,
    )

    # UV AOP
    m.fs.UV_aop = FlowsheetBlock(dynamic=False)
    build_uv_aop(m.fs.UV_aop, prop_package=m.fs.properties)

    # Decarbonator
    m.fs.decarbonator = FlowsheetBlock(dynamic=False)
    build_decarbonator(m.fs.decarbonator, prop_package=m.fs.properties)

    # Post-Treatment chemical addition units - ZO Models
    m.fs.post_treat_chem_list = [
        "calcium_hydroxide",
        "caustic",
        "sodium_hypochlorite",
        "sodium_bisulfite",
    ]

    for chem_name in m.fs.post_treat_chem_list:
        if chem_name == "sodium_hypochlorite":
            m.fs.add_component(
                chem_name + "_addition_post", FlowsheetBlock(dynamic=False)
            )
            build_chem_addition(
                m.fs.find_component(chem_name + "_addition_post"),
                chem_name,
                m.fs.properties,
            )
        else:
            m.fs.add_component(chem_name + "_addition", FlowsheetBlock(dynamic=False))
            build_chem_addition(
                m.fs.find_component(chem_name + "_addition"), chem_name, m.fs.properties
            )
    # Combined chemical list for operating conditions, scaling, and costing(?)
    m.fs.chemical_list = list(m.fs.pre_treat_chem_list) + list(
        m.fs.post_treat_chem_list
    )

    # Products and Disposal
    m.fs.disposal_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=["from_uf_disposal", "from_pro_brine", "from_tsro_brine"],
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)
    m.fs.disposal = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.disposal)

    return m


def add_wrd_connections(m):
    # Connect pre-UF chemical chain: feed -> chem1 -> chem2 -> ... -> UF
    for i, chem_name in enumerate((m.fs.pre_treat_chem_list)):
        unit = m.fs.find_component(chem_name + "_addition")
        if i == 0:
            a = Arc(
                source=m.fs.feed.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component(f"feed_to_{chem_name}", a)
        else:
            # Connect each chemical to the next
            prev = m.fs.pre_treat_chem_list[i - 1]
            prev_unit = m.fs.find_component(prev + "_addition")
            a = Arc(
                source=prev_unit.product.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component(f"{prev}_to_{chem_name}", a)

    last_chem = m.fs.find_component(m.fs.pre_treat_chem_list[-1] + "_addition")
    m.fs.pre_chem_to_uf_system = Arc(
        source=last_chem.product.outlet, destination=m.fs.uf_feed_separator.inlet
    )
    m.fs.uf_system_to_pro = Arc(
        source=m.fs.uf_product_mixer.outlet, destination=m.fs.ro_feed_separator.inlet
    )
    m.fs.pro_to_tsro_header = Arc(
        source=m.fs.ro_brine_mixer.outlet, destination=m.fs.tsro_header.inlet
    )
    m.fs.tsro_header_to_tsro_separator = Arc(
        source=m.fs.tsro_header.outlet, destination=m.fs.tsro_feed_separator.inlet
    )
    m.fs.pro_to_ro_system_product_mixer = Arc(
        source=m.fs.ro_product_mixer.outlet,
        destination=m.fs.ro_system_product_mixer.from_pro_product,
    )
    m.fs.pro_to_disposal_mixer = Arc(
        source=m.fs.ro_brine_mixer.outlet,
        destination=m.fs.disposal_mixer.from_pro_brine,
    )
    for t in m.fs.tsro_trains:
        a = Arc(
            source=m.fs.tsro_feed_separator.find_component(f"to_tsro{t}"),
            destination=m.fs.tsro_train[t].feed.inlet,
        )
        m.fs.add_component(f"tsro_header_to_tsro{t}", a)
        a = Arc(
            source=m.fs.tsro_train[t].product.outlet,
            destination=m.fs.ro_system_product_mixer.find_component(
                f"tsro{t}_to_ro_product"
            ),
        )
        m.fs.add_component(f"tsro{t}_to_ro_product", a)
        # m.fs.add_component(
        #     f"tsro{t}_to_product",
        #     Arc(
        #         source=m.fs.tsro_train[t].permeate.outlet,
        #         destination=m.fs.tsro_product_mixer.find_component("tsro" + str(t) + "_to_product").inlet,
        #     ),
        # )
        a = Arc(
            source=m.fs.tsro_train[t].disposal.outlet,
            destination=m.fs.tsro_brine_mixer.find_component(f"tsro{t}_to_brine"),
        )
        m.fs.add_component(f"tsro{t}_to_brine", a)

    m.fs.tsro_brine_mixer_to_disposal = Arc(
        source=m.fs.tsro_brine_mixer.outlet,
        destination=m.fs.disposal_mixer.from_tsro_brine,
    )
    m.fs.uf_disposal_to_disposal_mixer = Arc(
        source=m.fs.uf_disposal_mixer.outlet,
        destination=m.fs.disposal_mixer.from_uf_disposal,
    )
    m.fs.ro_system_product_mixer_to_uv = Arc(
        source=m.fs.ro_system_product_mixer.outlet, destination=m.fs.UV_aop.feed.inlet
    )
    m.fs.uv_to_decarbonator = Arc(
        source=m.fs.UV_aop.product.outlet, destination=m.fs.decarbonator.feed.inlet
    )

    # Chain post-treatment chemicals (decarbonator -> post1 -> post2 -> ... -> product)
    for i, chem_name in enumerate((m.fs.post_treat_chem_list)):
        unit = m.fs.find_component(chem_name + "_addition")
        if i == 0:
            # Connect decarb to first chemical
            a = Arc(
                source=m.fs.decarbonator.product.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component("decarb_to_" + chem_name, a)
        else:
            # Connect each chemical to the next
            prev = m.fs.post_treat_chem_list[i - 1]
            prev_unit = m.fs.find_component(prev + "_addition")
            a = Arc(
                source=prev_unit.product.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component(f"{prev}_to_{chem_name}", a)

    last_chem = m.fs.find_component(m.fs.post_treat_chem_list[-1] + "_addition")
    m.fs.post_chem_to_product = Arc(
        source=last_chem.product.outlet, destination=m.fs.product.inlet
    )
    m.fs.disposal_mixer_to_disposal = Arc(
        source=m.fs.disposal_mixer.outlet, destination=m.fs.disposal.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_inlet_conditions(m, Qin=2637*4, Cin=0.5, file="wrd_ro_inputs_8_19_21.yaml"):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + 27,
        },
        hold_state=True,
    )


def set_wrd_operating_conditions(m):

    # Operating conditions
    for chem_name in m.fs.chemical_list:
        set_chem_addition_op_conditions(
            blk=m.fs.find_component(chem_name + "_addition"), dose=0.1
        )

    set_uf_system_op_conditions(m)
    set_ro_system_op_conditions(m)
    
    for t in m.fs.tsro_trains:
        if t == m.fs.tsro_trains.first():
            m.fs.tsro_feed_separator.split_fraction[0, f"to_tsro{t}", "H2O"].set_value(1 / len(m.fs.tsro_trains))
        else:
            m.fs.tsro_feed_separator.split_fraction[0, f"to_tsro{t}", "H2O"].fix(1 / len(m.fs.tsro_trains))
        set_ro_stage_op_conditions(m.fs.tsro_train[t])

    set_uv_aop_op_conditions(m.fs.UV_aop)
    set_decarbonator_op_conditions(m.fs.decarbonator)


def set_wrd_system_scaling(m):
    # Properties Scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    # Does ZO property block also require scaling?
    # for chem_name in m.fs.chemical_list:
    #     set_chem_addition_scaling(blk=m.fs.find_component(chem_name + "_addition"))

    set_uf_system_scaling(m)
    set_ro_system_scaling(m)
    for t in m.fs.tsro_trains:
        set_ro_stage_scaling(m.fs.tsro_train[t])
    # add_UF_pump_scaling(m.fs.UF_pumps)
    # add_pump_scaling(m.fs.UF_pumps)
    # # add_separator_scaling(m.fs.UF) # Nothing to scale?
    # add_ro_scaling(m.fs.ro_system)
    add_uv_aop_scaling(m.fs.UV_aop)
    add_decarbonator_scaling(m.fs.decarbonator)


def initialize_wrd_system(m):
    m.fs.feed.initialize()
    # Initialize pre-UF chemical chain
    for i, chem_name in enumerate((m.fs.pre_treat_chem_list)):
        unit = m.fs.find_component(chem_name + "_addition")
        if i == 0:
            a = m.fs.find_component(f"feed_to_{chem_name}")
            propagate_state(a)
            initialize_chem_addition(unit)
        else:
            a = m.fs.find_component(f"{prev}_to_{chem_name}")
            propagate_state(a)
            initialize_chem_addition(unit)
        prev = chem_name
    last_chem = m.fs.find_component(m.fs.pre_treat_chem_list[-1] + "_addition")
    propagate_state(m.fs.pre_chem_to_uf_system)
    initialize_uf_system(m)
    propagate_state(m.fs.uf_system_to_pro)
    initialize_ro_system(
        m)
    propagate_state(m.fs.pro_to_tsro_header)
    m.fs.tsro_header.initialize()
    propagate_state(m.fs.tsro_header_to_tsro_separator)
    m.fs.tsro_feed_separator.initialize()
    # for t in m.fs.tsro_trains:
    #     a = m.fs.find_component(f"tsro_header_to_tsro{t}")
    #     propagate_state(a)
    #     initialize_ro_stage(m.fs.tsro_train[t])
    #     a = m.fs.ro_system_product_mixer.find_component(f"tsro{t}_to_ro_product")
    #     propagate_state(a)
    #     a = m.fs.tsro_brine_mixer.find_component(f"tsro{t}_to_brine")
    #     propagate_state(a)
        
        # a = Arc(
        #     source=m.fs.tsro_feed_separator.find_component(f"to_tsro{t}"),
        #     destination=m.fs.tsro_train[t].feed.inlet,
        # )
        # m.fs.add_component(f"tsro_header_to_tsro{t}", a)
        # a = Arc(
        #     source=m.fs.tsro_train[t].product.outlet,
        #     destination=m.fs.ro_system_product_mixer.find_component(
        #         f"tsro{t}_to_ro_product"
        #     ),
        # )
        # m.fs.add_component(f"tsro{t}_to_ro_product", a)
        # # m.fs.add_component(
        # #     f"tsro{t}_to_product",
        # #     Arc(
        # #         source=m.fs.tsro_train[t].permeate.outlet,
        # #         destination=m.fs.tsro_product_mixer.find_component("tsro" + str(t) + "_to_product").inlet,
        # #     ),
        # # )
        # a = Arc(
        #     source=m.fs.tsro_train[t].disposal.outlet,
        #     destination=m.fs.tsro_brine_mixer.find_component(f"tsro{t}_to_brine"),
        # )
        # m.fs.add_component(f"tsro{t}_to_brine", a)
    # propagate from last pre-UF chemical to UF
    # propagate_state(
    #     m.fs.find_component(m.fs.pre_treat_chem_list[-1] + "_to_translator")
    # )

    # propagate_state(m.fs.translator_to_uf_pumps)
    # # UF Pumps and separator
    # init_UF_pumps(m.fs.UF_pumps)
    # propagate_state(m.fs.uf_pumps_to_uf)
    # init_separator(m.fs.UF)

    # propagate_state(m.fs.UF_to_ro)
    # initialize_ro_system(m.fs.ro_system)
    # propagate_state(m.fs.ro_to_uv)
    # initialize_uv_aop(m.fs.UV_aop)
    # propagate_state(m.fs.uv_to_decarbonator)
    # initialize_decarbonator(m.fs.decarbonator)

    # # Initialize post-treatment chemical chain (downstream of decarbonator)
    # for i, chem_name in enumerate(m.fs.post_treat_chem_list):
    #     if i == 0:
    #         propagate_state(m.fs.find_component("decarb_to_" + chem_name))
    #     else:
    #         prev = m.fs.post_treat_chem_list[i - 1]
    #         propagate_state(m.fs.find_component(prev + "_to_" + chem_name))
    #     initialize_chem_addition(m.fs.find_component(chem_name + "_addition"))

    # propagate_state(m.fs.find_component(m.fs.post_treat_chem_list[-1] + "_to_product"))
    # m.fs.product.initialize()
    # propagate_state(m.fs.ro_waste_to_brine)
    # m.fs.brine.initialize()



def report_sj(sj, w=25):
    
    title = sj.name.replace("fs.", "").replace("_", " ").upper()

    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    flow_in = value(
        pyunits.convert(
            sj.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_in = value(
        pyunits.convert(
            sj.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f'{"INLET/OUTLET Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"INLET/OUTLET NaCl":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
    # flow_out = value(
    #     pyunits.convert(
    #         sj.outlet[0].flow_vol_phase["Liq"],
    #         to_units=pyunits.gallons / pyunits.minute,
    #     )
    # )
    # conc_out = value(
    #     pyunits.convert(
    #         sj.outlet[0].conc_mass_phase_comp["Liq", "NaCl"],
    #         to_units=pyunits.mg / pyunits.L,
    #     )
    # )
    # print(f'{"OUTLET Flow":<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}')
    # print(f'{"OUTLET NaCl":<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}')

def report_mixer(mixer, w=25):
    
    ms = mixer.find_component("mixed_state")
    title = mixer.name.replace("fs.", "").replace("_", " ").upper()
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    tot_flow_in = sum(
        value(
            value(
                pyunits.convert(
                    mixer.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                    to_units=pyunits.gallons / pyunits.minute,
                )
            )
        )
        for x in mixer.config.inlet_list
    )
    print(
        f'{"TOTAL INLET FLOW":<{w}s}{f"{tot_flow_in:<{w},.1f}"}{"gpm":<{w}s}'
    )
    for x in mixer.config.inlet_list:
        sb = mixer.find_component(f"{x}_state")
        flow_in = value(
            pyunits.convert(
                sb[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallons / pyunits.minute,
            )
        )
        tot_flow_in += flow_in
        conc_in = value(
            pyunits.convert(
                sb[0].conc_mass_phase_comp["Liq", "NaCl"],
                to_units=pyunits.mg / pyunits.L,
            )
        )
        print(
            f'{"   Flow " + x.replace("_", " ").title():<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}'
        )
        print(
            f'{"   NaCl " + x.replace("_", " ").title():<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}'
        )
    flow_out = value(
        pyunits.convert(
            ms[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_out = value(
        pyunits.convert(
            ms[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f'{"Outlet Flow":<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"Outlet NaCl":<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}')

def report_separator(sep, w=25):

    title = sep.name.replace("fs.", "").replace("_", " ").upper()

    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    ms = sep.find_component("mixed_state")
    flow_in = value(
        pyunits.convert(
            ms[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_in = value(
        pyunits.convert(
            ms[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f'{"INLET Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"INLET NaCl":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
    tot_flow_out = sum(
        value(
            value(
                pyunits.convert(
                    sep.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                    to_units=pyunits.gallons / pyunits.minute,
                )
            )
        )
        for x in sep.config.outlet_list
    )
    print(f'{"TOTAL OUTLET FLOW":<{w}s}{f"{tot_flow_out:<{w},.1f}"}{"gpm":<{w}s}')
    for x in sep.config.outlet_list:
        sb = sep.find_component(f"{x}_state")
        flow_out = value(
            pyunits.convert(
                sb[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallons / pyunits.minute,
            )
        )
        conc_out = value(
            pyunits.convert(
                sb[0].conc_mass_phase_comp["Liq", "NaCl"],
                to_units=pyunits.mg / pyunits.L,
            )
        )
        print(
            f'{"   Flow " + x.replace("_", " ").title():<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}'
        )
        print(
            f'{"   NaCl " + x.replace("_", " ").title():<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}'
        )

def main(num_pro_trains=4, num_tsro_trains=None, num_stages=2):
    m = build_wrd_system(
        num_pro_trains=num_pro_trains,
        num_tsro_trains=num_tsro_trains,
        num_stages=num_stages,
    )
    # assert_units_consistent(m)
    add_wrd_connections(m)
    # print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    # set_inlet_conditions(m)
    set_wrd_operating_conditions(m)
    set_wrd_system_scaling(m)
    calculate_scaling_factors(m)
    set_inlet_conditions(m)
    initialize_wrd_system(m)

    # print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    print(f"dof = {degrees_of_freedom(m)}")

    # initialize_wrd_system(m)
    # try:
    #     results = solve(m)
    #     assert_optimal_termination(results)
    # except:
    #     print_infeasible_constraints(m)
    #     print("\n--------- Failed to Solve ---------\n")
    return m

if __name__ == "__main__":
    num_stages = 2
    
    m = main(num_stages=num_stages)
    # from wrd.components.UF_separator import report_uf_separator
    # m.fs.tsro_header.properties[0].display()
    # report_ro_system(m)
    report_mixer(m.fs.uf_product_mixer)
    report_separator(m.fs.ro_feed_separator)
    report_mixer(m.fs.ro_product_mixer)
    report_mixer(m.fs.ro_brine_mixer)
    report_sj(m.fs.tsro_header)
    # report_ro_train(m.fs.train[1])
    report_separator(m.fs.tsro_feed_separator)

    # m = build_wrd_system(num_stages=num_stages, date=date)
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
