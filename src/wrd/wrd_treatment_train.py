from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
    value,
    TransformationFactory,
)

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.models.unit_models import Product, Feed

from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.costing import WaterTAPCosting

from wrd.components.decarbonator import *
from wrd.components.uv_aop import *
from wrd.components.pump import *
from wrd.components.UF_system import *
from wrd.components.ro_system import *
from wrd.components.ro_stage import *
from wrd.components.ro import report_ro
from wrd.components.chemical_addition import *
from wrd.components.brine_disposal import *
from wrd.utilities import *
from srp.utils import touch_flow_and_conc
from models import HeadLoss, Source


def build_wrd_system(num_pro_trains=4, num_tsro_trains=None, num_stages=2, file=None):

    if file is None:
        raise ValueError("Input file must be provided to build WRD system.")

    if num_tsro_trains is None:
        num_tsro_trains = num_pro_trains

    m = ConcreteModel()
    m.num_pro_trains = num_pro_trains
    m.num_tsro_trains = num_tsro_trains
    m.num_stages = num_stages
    m.fs = FlowsheetBlock(dynamic=False)

    config = get_config_file(file)
    m.fs.config_data = load_config(config)

    config = get_config_file("chemical_addition.yaml")
    m.fs.chem_data = load_config(config)

    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # configure costing parameters
    m.fs.costing.base_currency = pyunits.USD_2021
    m.fs.costing.base_period = pyunits.year
    m.fs.costing.utilization_factor.fix(1)
    m.fs.costing.maintenance_labor_chemical_factor.fix(0)
    m.fs.costing.electricity_cost.fix(0.15)

    # Add units
    m.fs.feed = Source(property_package=m.fs.properties)
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
    build_uf_system(
        m=m, num_trains=num_pro_trains, prop_package=m.fs.properties, file=file
    )

    # PRO System
    build_ro_system(
        m=m,
        num_trains=num_pro_trains,
        num_stages=num_stages,
        prop_package=m.fs.properties,
        file=file,
    )

    # TSRO System
    m.fs.tsro_trains = Set(initialize=range(1, num_tsro_trains + 1))
    # m.fs.tsro_header = StateJunction(property_package=m.fs.properties)
    m.fs.tsro_header = HeadLoss(property_package=m.fs.properties)
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
        build_ro_stage(
            m.fs.tsro_train[t], stage_num=3, prop_package=m.fs.properties, file=file
        )

    ro_system_prod_mixer_inlet_list = ["from_pro_product"] + [
        f"tsro{t}_to_ro_product" for t in m.fs.tsro_trains
    ]

    m.fs.ro_system_product_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=ro_system_prod_mixer_inlet_list,
        momentum_mixing_type=MomentumMixingType.none,
    )
    touch_flow_and_conc(m.fs.ro_system_product_mixer)

    m.fs.tsro_brine_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=[f"tsro{t}_to_brine" for t in m.fs.tsro_trains],
        momentum_mixing_type=MomentumMixingType.none,
    )
    touch_flow_and_conc(m.fs.tsro_brine_mixer)

    m.fs.total_tsro_pump_power = Expression(
        expr=sum(
            m.fs.tsro_train[i].pump.unit.work_mechanical[0] for i in m.fs.tsro_trains
        )
    )

    # UV AOP
    m.fs.UV_aop = FlowsheetBlock(dynamic=False)
    build_uv_aop(m.fs.UV_aop, prop_package=m.fs.properties)

    # Decarbonator
    m.fs.decarbonator = FlowsheetBlock(dynamic=False)
    build_decarbonator(m.fs.decarbonator, prop_package=m.fs.properties)

    # Post-Treatment chemical addition units
    m.fs.post_treat_chem_list = [
        "calcium_hydroxide",
        "sodium_hydroxide",
        "sodium_bisulfite",
    ]

    for chem_name in m.fs.post_treat_chem_list:
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
    touch_flow_and_conc(m.fs.disposal_mixer)

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)

    # Brine disposal
    m.fs.disposal = FlowsheetBlock(dynamic=False)
    build_brine_disposal(m.fs.disposal, prop_package=m.fs.properties, file=file)

    # Overall System Metrics
    m.fs.system_recovery = Expression(
        expr=m.fs.product.properties[0].flow_vol_phase["Liq"]
        / m.fs.feed.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.total_system_pump_power = Expression(
        expr=pyunits.convert(
            m.fs.total_uf_pump_power
            + m.fs.total_ro_pump_power
            + m.fs.total_tsro_pump_power,
            to_units=pyunits.kW,
        )
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
        m.fs.add_component(f"tsro_separator_to_tsro{t}", a)
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

    m.fs.ro_system_product_mixer_to_uv = Arc(
        source=m.fs.ro_system_product_mixer.outlet, destination=m.fs.UV_aop.feed.inlet
    )
    m.fs.uv_to_decarbonator = Arc(
        source=m.fs.UV_aop.product.outlet, destination=m.fs.decarbonator.feed.inlet
    )

    # Chain post-treatment chemicals (decarbonator -> post1 -> post2 -> ... -> product)
    for i, chem_name in enumerate(m.fs.post_treat_chem_list):
        if chem_name == "sodium_hypochlorite":
            unit = m.fs.find_component(chem_name + "_addition_post")
        else:
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
        source=m.fs.disposal_mixer.outlet, destination=m.fs.disposal.feed.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_wrd_inlet_conditions(m, Qin=None, Cin=None, Tin=None):
    # IMO it makes sense Qin to be from the yaml for the wrd flowsheet so all case study inputs are in one place.
    # Other component files Q and C can be hard coded for testing.
    if Qin is None:
        Qin = get_config_value(m.fs.config_data, "feed_flow_water", "feed_stream")
    else:
        Qin = Qin * pyunits.gallons / pyunits.minute

    if Cin is None:
        Cin = get_config_value(
            m.fs.config_data, "feed_conductivity", "feed_stream"
        ) * get_config_value(
            m.fs.config_data, "feed_conductivity_conversion", "feed_stream"
        )
    else:
        Cin = Cin * pyunits.g / pyunits.L

    if Tin is None:
        Tin = get_config_value(m.fs.config_data, "feed_temperature", "feed_stream")
    else:
        Tin = Tin * pyunits.K

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): (Qin * m.num_pro_trains),
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin,
            ("temperature", None): Tin,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )


def set_wrd_operating_conditions(m):

    # Operating conditions
    for chem_name in m.fs.chemical_list:
        set_chem_addition_op_conditions(
            blk=m.fs.find_component(chem_name + "_addition")
        )

    set_uf_system_op_conditions(m)

    set_ro_system_op_conditions(m)

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

    set_uv_aop_op_conditions(m.fs.UV_aop)

    set_decarbonator_op_conditions(m.fs.decarbonator)
    m.fs.tsro_header.TSRO_header_loss = get_config_value(
        m.fs.config_data, "header_loss", "reverse_osmosis_1d", "stage_3"
    )
    m.fs.tsro_header.control_volume.deltaP[0].fix(m.fs.tsro_header.TSRO_header_loss)
    m.fs.ro_system_product_mixer.outlet.pressure[0].fix(101325)
    m.fs.tsro_brine_mixer.outlet.pressure[0].fix(101325)
    m.fs.disposal_mixer.outlet.pressure[0].fix(101325)


def set_wrd_system_scaling(m):

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

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

    m.fs.uf_disposal_mixer.initialize()
    propagate_state(m.fs.uf_disposal_to_disposal_mixer)

    propagate_state(m.fs.uf_system_to_pro)
    initialize_ro_system(m)

    propagate_state(m.fs.pro_to_ro_system_product_mixer)

    propagate_state(m.fs.pro_to_tsro_header)

    m.fs.tsro_header.initialize()
    propagate_state(m.fs.tsro_header_to_tsro_separator)
    m.fs.tsro_feed_separator.initialize()

    for t in m.fs.tsro_trains:
        a = m.fs.find_component(f"tsro_separator_to_tsro{t}")
        propagate_state(a)
        initialize_ro_stage(m.fs.tsro_train[t])
        a = m.fs.find_component(f"tsro{t}_to_ro_product")
        propagate_state(a)
        a = m.fs.find_component(f"tsro{t}_to_brine")
        propagate_state(a)

    m.fs.tsro_brine_mixer.initialize()
    propagate_state(m.fs.tsro_brine_mixer_to_disposal)

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
            propagate_state(a)
            initialize_chem_addition(unit)
        else:
            a = m.fs.find_component(f"{prev}_to_{chem_name}")
            propagate_state(a)
            initialize_chem_addition(unit)
        prev = chem_name
    propagate_state(m.fs.post_chem_to_product)
    m.fs.product.initialize()

    m.fs.disposal_mixer.initialize()
    propagate_state(m.fs.disposal_mixer_to_disposal)

    initialize_brine_disposal(m.fs.disposal)


def add_wrd_system_costing(m, source_cost=0.15):

    m.fs.feed.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.source.unit_cost.fix(source_cost)

    add_uf_system_costing(m, costing_package=m.fs.costing)
    add_ro_system_costing(m, costing_package=m.fs.costing)
    cost_uv_aop(m.fs.UV_aop, costing_package=m.fs.costing)
    add_brine_disposal_costing(m.fs.disposal, costing_package=m.fs.costing)
    cost_decarbonator(m.fs.decarbonator, costing_package=m.fs.costing)
    for chem_name in m.fs.chemical_list:
        add_chem_addition_costing(
            blk=m.fs.find_component(chem_name + "_addition"),
            costing_package=m.fs.costing,
        )

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.properties[0].flow_vol_phase["Liq"],
        name="SEC",
    )
    m.fs.costing.initialize()


def report_tsro(m, w=30):
    title = "TSRO Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    for t in m.fs.tsro_trains:
        tsro_stage = m.fs.tsro_train[t]
        report_pump(tsro_stage.pump)
        report_ro(tsro_stage.ro, w=w)


def report_wrd_comparison_metrics(m, w=30):

    print("Comparative Metrics:")
    # Print UF system metrics
    title = f" UF System Metrics"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "-" * side + f" {title} " + "-" * side
    print(f"\n{header}\n")
    for i in m.fs.uf_trains:
        print(
            f'{f"UF Train {i} Feed Flow":<{w}s}{value(pyunits.convert(m.fs.uf_train[i].feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"UF Train {i} Pump Power":<{w}s}{value(pyunits.convert(m.fs.uf_train[i].pump.unit.work_mechanical[0], to_units=pyunits.kW)):<{w}.3f}{"kW"}'
        )

    # Print PRO train metrics
    for i in m.fs.trains:
        title = f"PRO Train {i} Metrics"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "-" * side + f" {title} " + "-" * side
        print(f"\n{header}\n")

        # Print stage-by-stage pump powers and permeate flows
        for j in m.fs.train[i].stages:
            print(
                f'{f"  Stage {j} Pump Power":<{w}s}{value(pyunits.convert(m.fs.train[i].stage[j].pump.unit.work_mechanical[0], to_units=pyunits.kW)):<{w}.3f}{"kW"}'
            )
            print(
                f'{f"  Stage {j} Perm Flow":<{w}s}{value(pyunits.convert(m.fs.train[i].stage[j].product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
            )
            print(
                f'{f"Perm Conc":<{w}s}{value(pyunits.convert(m.fs.train[i].stage[j].ro.unit.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.liter)):<{w}.3f}{f"mg/L"}'
            )
            print(
                f'{f"  Stage {j} Recovery":<{w}s}{value(m.fs.train[i].stage[j].ro.unit.recovery_vol_phase[0, "Liq"])*100:<{w}.3f}{"%"}'
            )

        # Print train totals
        print(
            f'{f"  Train {i} Total Pump Power":<{w}s}{value(pyunits.convert(m.fs.train[i].total_pump_power, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
        )
        print(
            f'{f"  Train {i} Total Perm Flow":<{w}s}{value(pyunits.convert(m.fs.train[i].product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )

    # Print TSRO train metrics
    for t in m.fs.tsro_trains:
        title = f"TSRO Train {t} Metrics"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "-" * side + f" {title} " + "-" * side
        print(f"\n{header}\n")

        print(
            f'{f"  TSRO {t} Pump Power":<{w}s}{value(pyunits.convert(m.fs.tsro_train[t].pump.unit.work_mechanical[0], to_units=pyunits.kW)):<{w}.3f}{"kW"}'
        )
        print(
            f'{f"  TSRO {t} Perm Flow":<{w}s}{value(pyunits.convert(m.fs.tsro_train[t].product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"  TSRO {t} Recovery":<{w}s}{value(m.fs.tsro_train[t].ro.unit.recovery_vol_phase[0, "Liq"])*100:<{w}.3f}{"%"}'
        )
    # Print decarb and UV metrics
    title = f"Decarbonator and UV AOP Metrics"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "-" * side + f" {title} " + "-" * side
    print(f"\n{header}\n")

    print(
        f'{f"UV AOP Feed Flow":<{w}s}{value(pyunits.convert(m.fs.UV_aop.feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"UV AOP Energy Use":<{w}s}{value(pyunits.convert(m.fs.UV_aop.unit.power_consumption, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    print(
        f'{f"Decarbonator Energy Use":<{w}s}{value(pyunits.convert(m.fs.decarbonator.unit.power_consumption, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )

    # Costs
    if m.fs.find_component("costing") is not None:
        title = f"Flow Costs"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "-" * side + f" {title} " + "-" * side
        print(f"\n{header}\n")
        # print(
        #     f'{f"Levelized Cost of Water":<{w}s}{value(pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2021  / pyunits.m**3)):<{w}.3f}{"$/m3"}'
        # )
        for key in m.fs.costing.aggregate_flow_costs:
            print(
                f'{f"{key}":<{w}s}{value(m.fs.costing.aggregate_flow_costs[key]):<{w}.3f}{"$/yr"}'
            )
        print(
            f'{f"Brine Disposal Opex":<{w}s}{value(m.fs.disposal.unit.costing.variable_operating_cost):<{w}.2f}{"$/yr"}'
        )
        print(
            f'{f"Feed Opex":<{w}s}{value(m.fs.feed.costing.variable_operating_cost):<{w}.3f}{"$/yr"}'
        )


def report_wrd(m, w=30, add_comp_metrics=False):

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
        m.fs.disposal.feed.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gallon / pyunits.minute,
    )
    brine_flow_mgd = pyunits.convert(
        brine_flow, to_units=pyunits.Mgallons / pyunits.day
    )
    brine_conc = pyunits.convert(
        m.fs.disposal.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        to_units=pyunits.mg / pyunits.L,
    )

    title = f"WRD System Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "/" * side + f" {title} " + "\\" * side
    sep = "\n" + "@" * (3 * w) + "\n"
    print(f"\n\n{header}\n")

    for i, chem_name in enumerate(m.fs.pre_treat_chem_list):
        unit = m.fs.find_component(chem_name + "_addition")
        report_chem_addition(unit, w=w)

    print(sep)
    report_uf_system(m, w=w)
    print(sep)
    report_ro_system(m, w=w)
    print(sep)

    report_mixer(m.fs.ro_brine_mixer, w=w)
    print(sep)
    report_head_loss(m.fs.tsro_header, w=w)
    print(sep)

    for t in m.fs.tsro_trains:
        title = f"TSRO Stage Report - Train {t}"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "=" * side + f" {title} " + "=" * side
        print(f"\n{header}\n")
        report_ro_stage(m.fs.tsro_train[t], w=w)

    print(sep)
    report_uv(m.fs.UV_aop, w=w)
    print(sep)
    report_decarbonator(m.fs.decarbonator, w=w)
    print(sep)

    for chem_name in m.fs.post_treat_chem_list:
        unit = m.fs.find_component(chem_name + "_addition")
        report_chem_addition(unit, w=w)

    print(sep)
    report_mixer(m.fs.tsro_brine_mixer, w=w)
    report_mixer(m.fs.uf_disposal_mixer, w=w)
    print(sep)
    report_mixer(m.fs.ro_system_product_mixer, w=w)
    report_mixer(m.fs.disposal_mixer, w=w)
    print(sep)
    report_brine_disposal(m.fs.disposal, w=w)
    print(sep)

    title = f"WRD System Summary"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "/" * side + f" {title} " + "\\" * side
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
    print()
    print(
        f'{f"Total Pumping Power":<{w}s}{value(pyunits.convert(m.fs.total_system_pump_power, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )
    print(sep)
    if add_comp_metrics:
        report_wrd_comparison_metrics(m, w=w)


def main(
    num_pro_trains=4,
    num_tsro_trains=None,
    num_pro_stages=2,
    file=None,
):

    m = build_wrd_system(
        num_pro_trains=num_pro_trains,
        num_tsro_trains=num_tsro_trains,
        num_stages=num_pro_stages,
        file=file,
    )
    add_wrd_connections(m)
    set_wrd_system_scaling(m)
    calculate_scaling_factors(m)
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_wrd_inlet_conditions(m)
    set_wrd_operating_conditions(m)
    print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    assert degrees_of_freedom(m) == 0
    initialize_wrd_system(m)

    add_wrd_system_costing(m)

    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_wrd(m, add_comp_metrics=True)

    return m


if __name__ == "__main__":
    num_pro_trains = 1
    file = "wrd_inputs_3_13_21.yaml"
    m = main(num_pro_trains=num_pro_trains, file=file)
