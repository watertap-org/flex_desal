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
from wrd.components.ro import report_ro
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

    # Connections by pass this mixer, so can be removed
    # m.fs.tsro_product_mixer = Mixer(
    #     property_package=m.fs.properties,
    #     inlet_list=[f"tsro{i}_to_product" for i in m.fs.tsro_trains],
    #     momentum_mixing_type=MomentumMixingType.none,
    # )

    ro_system_prod_mixer_inlet_list = ["from_pro_product"] + [
        f"tsro{t}_to_ro_product" for t in m.fs.tsro_trains
    ]

    m.fs.ro_system_product_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=ro_system_prod_mixer_inlet_list,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.tsro_brine_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=[f"tsro{t}_to_brine" for t in m.fs.tsro_trains],
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
        "sodium_hydroxide",
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
    m.fs.chemical_list = m.fs.pre_treat_chem_list + m.fs.post_treat_chem_list

    # Products and Disposal
    m.fs.disposal_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=["from_uf_disposal", "from_tsro_brine"],
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)
    m.fs.disposal = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.disposal)

    m.fs.system_recovery = Expression(
        expr=m.fs.product.properties[0].flow_vol_phase["Liq"]
        / m.fs.feed.properties[0].flow_vol_phase["Liq"]
    )

    return m


def add_wrd_connections(m):
    # Connect pre-UF chemical chain: feed -> chem1 -> chem2 -> ... -> UF
    for i, chem_name in enumerate(m.fs.pre_treat_chem_list):
        unit = m.fs.find_component(chem_name + "_addition")
        if i == 0:
            a = Arc(
                source=m.fs.feed.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component(f"feed_to_{chem_name}", a)
        else:
            # Connect each chemical to the next
            prev_unit = m.fs.find_component(f"{prev}_addition")
            a = Arc(
                source=prev_unit.product.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component(f"{prev}_to_{chem_name}", a)
        prev = chem_name

    last_chem = m.fs.find_component(f"{prev}_addition")
    m.fs.pre_chem_to_uf_system = Arc(
        source=last_chem.product.outlet, destination=m.fs.uf_feed_separator.inlet
    )
    m.fs.uf_system_to_pro = Arc(
        source=m.fs.uf_product_mixer.outlet, destination=m.fs.ro_feed_separator.inlet
    )
    #####
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
    # There is no direct line from pro to brine
    # m.fs.pro_to_disposal_mixer = Arc(
    #     source=m.fs.ro_brine_mixer.outlet,
    #     destination=m.fs.disposal_mixer.from_pro_brine,
    # )
    m.fs.ro_system_product_mixer_to_uv = Arc(
        source=m.fs.ro_system_product_mixer.outlet, destination=m.fs.UV_aop.feed.inlet
    )
    m.fs.uv_to_decarbonator = Arc(
        source=m.fs.UV_aop.product.outlet, destination=m.fs.decarbonator.feed.inlet
    )

    # Chain post-treatment chemicals (decarbonator -> post1 -> post2 -> ... -> product)
    for i, chem_name in enumerate(m.fs.post_treat_chem_list):
        unit = m.fs.find_component(f"{chem_name}_addition")
        if i == 0:
            # Connect decarb to first chemical
            a = Arc(
                source=m.fs.decarbonator.product.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component(f"decarb_to_{chem_name}", a)
        else:
            # Connect each chemical to the next
            prev_unit = m.fs.find_component(f"{prev}_addition")
            a = Arc(
                source=prev_unit.product.outlet,
                destination=unit.feed.inlet,
            )
            m.fs.add_component(f"{prev}_to_{chem_name}", a)
        prev = chem_name

    last_chem = m.fs.find_component(f"{prev}_addition")
    m.fs.post_chem_to_product = Arc(
        source=last_chem.product.outlet, destination=m.fs.product.inlet
    )
    m.fs.disposal_mixer_to_disposal = Arc(
        source=m.fs.disposal_mixer.outlet, destination=m.fs.disposal.inlet
    )
    ##### Is this TSRO Bypass?
    # m.fs.ro_to_disposal = Arc(
    #     source=m.fs.tsro_header.outlet, destination=m.fs.disposal.inlet
    # )

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
    print(degrees_of_freedom(m))


def set_wrd_operating_conditions(m):
    # Operating conditions
    for chem_name in m.fs.chemical_list:
        # Load dose from yaml??
        set_chem_addition_op_conditions(
            blk=m.fs.find_component(chem_name + "_addition"), dose=0.1
        )
        if chem_name == "sodium_hypochlorite":
            set_chem_addition_op_conditions(
            blk=m.fs.find_component(chem_name + "_addition_post"), dose=0.1
        )
    print(degrees_of_freedom(m))
    set_uf_system_op_conditions(m)
    print(degrees_of_freedom(m))
    set_ro_system_op_conditions(m)
    print(degrees_of_freedom(m))
    for t in m.fs.tsro_trains:
        if t != m.fs.tsro_trains.first():
            m.fs.tsro_feed_separator.split_fraction[0, f"to_tsro{t}", "H2O"].fix(
                1 / len(m.fs.tsro_trains)
            )
            m.fs.tsro_feed_separator.split_fraction[0, f"to_tsro{t}", "NaCl"].fix(
                1 / len(m.fs.tsro_trains)
            )
        else:
            m.fs.tsro_feed_separator.split_fraction[0, f"to_tsro{t}", "H2O"].set_value(
                1 / len(m.fs.tsro_trains)
            )
            m.fs.tsro_feed_separator.split_fraction[0, f"to_tsro{t}", "NaCl"].set_value(
                1 / len(m.fs.tsro_trains)
            )
        set_ro_stage_op_conditions(m.fs.tsro_train[t])
    print(degrees_of_freedom(m))
    set_uv_aop_op_conditions(m.fs.UV_aop)
    print(degrees_of_freedom(m))
    set_decarbonator_op_conditions(m.fs.decarbonator)
    print(degrees_of_freedom(m))

    # m.fs.tsro_product_mixer.outlet.pressure[0].fix(101325) #Removed mixer  
    m.fs.ro_system_product_mixer.outlet.pressure[0].fix(101325)
    m.fs.tsro_brine_mixer.outlet.pressure[0].fix(101325)
    m.fs.disposal_mixer.outlet.pressure[0].fix(101325)

    


def set_wrd_system_scaling(m):
    # Properties Scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    # NO Chem Addition Scaling?
    set_uf_system_scaling(m)
    set_ro_system_scaling(m)
    for t in m.fs.tsro_trains:
        set_ro_stage_scaling(m.fs.tsro_train[t])
    add_uv_aop_scaling(m.fs.UV_aop)
    add_decarbonator_scaling(m.fs.decarbonator)


def initialize_wrd_system(m):

    m.fs.feed.initialize()

    # Initialize pre-UF chemical chain
    for i, chem_name in enumerate(m.fs.pre_treat_chem_list):
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

    propagate_state(m.fs.pre_chem_to_uf_system)
    initialize_uf_system(m)
    propagate_state(m.fs.uf_system_to_pro)
    initialize_ro_system(m)

    # propagate_state(m.fs.ro_to_product)
    # m.fs.product.initialize()
    # propagate_state(m.fs.ro_to_disposal)
    
    propagate_state(m.fs.pro_to_ro_system_product_mixer)
    # propagate_state(m.fs.pro_to_disposal_mixer)
    propagate_state(m.fs.pro_to_tsro_header)
    m.fs.tsro_header.initialize()
    propagate_state(m.fs.tsro_header_to_tsro_separator)
    m.fs.tsro_feed_separator.initialize()
    for t in m.fs.tsro_trains:
        a = m.fs.find_component(f"tsro_header_to_tsro{t}")
        propagate_state(a)
        initialize_ro_stage(m.fs.tsro_train[t])
        a = m.fs.find_component(f"tsro{t}_to_ro_product")
        propagate_state(a)
        a = m.fs.find_component(f"tsro{t}_to_brine")
        propagate_state(a)

    m.fs.tsro_brine_mixer.initialize()
    propagate_state(m.fs.tsro_brine_mixer_to_disposal) # This mixer could be removed?

    m.fs.uf_disposal_mixer.initialize()
    propagate_state(m.fs.uf_disposal_to_disposal_mixer)

    m.fs.ro_system_product_mixer.initialize()
    propagate_state(m.fs.ro_system_product_mixer_to_uv)

    initialize_uv_aop(m.fs.UV_aop)
    propagate_state(m.fs.uv_to_decarbonator)

    initialize_decarbonator(m.fs.decarbonator)

    # Chain post-treatment chemicals (decarbonator -> post1 -> post2 -> ... -> product)
    for i, chem_name in enumerate((m.fs.post_treat_chem_list)):
        unit = m.fs.find_component(f"{chem_name}_addition")
        if i == 0:
            # Connect decarb to first chemical
            a = m.fs.find_component(f"decarb_to_{chem_name}")
            # m.fs.add_component("decarb_to_" + chem_name, a)
            propagate_state(a)
            initialize_chem_addition(unit)
        else:
            a = m.fs.find_component(f"{prev}_to_{chem_name}")
            propagate_state(a)
            initialize_chem_addition(unit)
        prev = chem_name
    propagate_state(m.fs.post_chem_to_product)
    # m.fs.post_chem_to_product = Arc(
    #     source=last_chem.product.outlet, destination=m.fs.product.inlet
    # )
    m.fs.product.initialize()

    m.fs.disposal_mixer.initialize()
    propagate_state(m.fs.disposal_mixer_to_disposal)

    m.fs.disposal.initialize()


def report_sj(sj, w=25):
    #sj = state junction
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
            pyunits.convert(
                mixer.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallons / pyunits.minute,
            )
        )
        for x in mixer.config.inlet_list
    )
    print(f'{"TOTAL INLET FLOW":<{w}s}{f"{tot_flow_in:<{w},.1f}"}{"gpm":<{w}s}')
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


def report_tsro(m, w=30):
    title = "TSRO Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    for t in m.fs.tsro_trains:
        tsro_stage = m.fs.m.fs.sro_train[t]
        report_pump(tsro_stage.pump)
        report_ro(tsro_stage.ro, w=w)

def report_wrd(m, w=30):

    feed_flow = pyunits.convert(
        m.fs.feed.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.minute,
    )
    feed_flow_mgd = pyunits.convert(feed_flow, to_units=pyunits.Mgallons / pyunits.day)
    feed_conc = pyunits.convert(
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        to_units=pyunits.mg / pyunits.L,
    )

    prod_flow = pyunits.convert(
        m.fs.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.minute,
    )
    prod_flow_mgd = pyunits.convert(prod_flow, to_units=pyunits.Mgallons / pyunits.day)
    prod_conc = pyunits.convert(
        m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        to_units=pyunits.mg / pyunits.L,
    )

    brine_flow = pyunits.convert(
        m.fs.disposal.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.minute,
    )
    brine_flow_mgd = pyunits.convert(
        brine_flow, to_units=pyunits.Mgallons / pyunits.day
    )
    brine_conc = pyunits.convert(
        m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        to_units=pyunits.mg / pyunits.L,
    )

    title = f"WRD Report Performance"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "_" * side + f" {title} " + "_" * side
    print(f"\n\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{"System Recovery":<{w}s}{value(m.fs.system_recovery)*100:<{w}.3f}{"%":<{w}s}'
    )
    print(f'{f"Feed Flow (gpm)":<{w}s}{value(feed_flow):<{w}.3f}{"gpm"}')
    print(f'{f"Feed Flow (MGD)":<{w}s}{value(feed_flow_mgd):<{w}.3f}{"MGD"}')
    print(f'{f"Feed Conc":<{w}s}{value(feed_conc):<{w}.3f}{"mg/L"}')
    print(f'{f"Product Flow (gpm)":<{w}s}{value(prod_flow):<{w}.3f}{"gpm"}')
    print(f'{f"Product Flow (MGD)":<{w}s}{value(prod_flow_mgd):<{w}.3f}{"MGD"}')
    print(f'{f"Product Conc":<{w}s}{value(prod_conc):<{w}.3f}{"mg/L"}')
    print(f'{f"Brine Flow (gpm)":<{w}s}{value(brine_flow):<{w}.3f}{"gpm"}')
    print(f'{f"Brine Flow (MGD)":<{w}s}{value(brine_flow_mgd):<{w}.3f}{"MGD"}')
    print(f'{f"Brine Conc":<{w}s}{value(brine_conc):<{w}.3f}{"mg/L"}')

    # for i, chem_name in enumerate(m.fs.pre_treat_chem_list):
    #     unit = m.fs.find_component(chem_name + "_addition")
    #     report_chem_addition(unit, w=w )

    # report_uf_system(m, w=w)
    report_ro_system(m, w=w)
    report_mixer(m.fs.ro_brine_mixer)
    report_sj(m.fs.tsro_header)
    report_tsro(m,w=w)


def main(num_pro_trains=4, num_tsro_trains=4, num_stages=2):
    m = build_wrd_system(
        num_pro_trains=num_pro_trains,
        num_tsro_trains=num_tsro_trains,
        num_stages=num_stages,
    )
    add_wrd_connections(m)
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_inlet_conditions(m,Qin=2637*num_pro_trains)
    set_wrd_operating_conditions(m)
    print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    set_wrd_system_scaling(m)
    calculate_scaling_factors(m)
    assert degrees_of_freedom(m) == 0
    initialize_wrd_system(m)
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
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_wrd(m)
    # from wrd.components.UF_separator import report_uf_separator
    # m.fs.tsro_header.properties[0].display()
    # report_ro_system(m)
    # report_mixer(m.fs.uf_product_mixer)
    # report_separator(m.fs.ro_feed_separator)
    # report_mixer(m.fs.ro_product_mixer)
    # report_mixer(m.fs.ro_brine_mixer)
    # report_sj(m.fs.tsro_header)
    # # report_ro_train(m.fs.train[1])
    # report_separator(m.fs.tsro_feed_separator)
    # x = pyunits.convert(m.fs.feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)
    # print(f"\nFeed Flowrate: {x():.1f} gpm\n")

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
