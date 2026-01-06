from pyomo.environ import (
    ConcreteModel,
    Expression,
    Set,
    value,
    assert_optimal_termination,
    units as pyunits,
    value,
    TransformationFactory,
)
from pyomo.network import Arc

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    MomentumMixingType,
    Mixer,
    Separator,
    Feed,
    Product,
    SplittingType,
)
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.costing import WaterTAPCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from wrd.components.ro_train import *
from wrd.components.pump import report_pump
from srp.utils import touch_flow_and_conc

solver = get_solver()


def build_ro_system(
    m=None,
    num_trains=3,
    num_stages=3,
    split_fractions=None,
    prop_package=None,
    file="wrd_inputs_8_19_21.yaml",
):

    if m is None:
        m = ConcreteModel()
        m.standalone = True
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = NaClParameterBlock()
        m.fs.feed = Feed(property_package=m.fs.properties)
        touch_flow_and_conc(m.fs.feed)
        m.fs.product = Product(property_package=m.fs.properties)
        touch_flow_and_conc(m.fs.product)
        m.fs.disposal = Product(property_package=m.fs.properties)
        touch_flow_and_conc(m.fs.disposal)

    else:
        m.standalone = False

    if prop_package is None:
        prop_package = m.fs.properties

    m.num_trains = num_trains
    m.num_stages = num_stages
    m.fs.trains = Set(initialize=range(1, m.num_trains + 1))
    m.fs.train = FlowsheetBlock(m.fs.trains, dynamic=False)

    outlet_list = [f"train{i}" for i in m.fs.trains]

    m.fs.ro_feed_separator = Separator(
        property_package=m.fs.properties,
        outlet_list=outlet_list,
        split_basis=SplittingType.componentFlow,
    )

    if split_fractions is None:
        m.fs.ro_feed_separator.even_split = 1.0 / len(outlet_list)
    else:
        raise NotImplementedError("Custom split fractions not yet implemented.")

    perm_inlet_list = [f"perm_inlet{i}" for i in m.fs.trains]

    m.fs.ro_product_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.minimize,
        inlet_list=perm_inlet_list,
    )
    touch_flow_and_conc(m.fs.ro_product_mixer)

    brine_inlet_list = [f"brine_inlet{i}" for i in m.fs.trains]

    m.fs.ro_brine_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.minimize,
        inlet_list=brine_inlet_list,
    )
    touch_flow_and_conc(m.fs.ro_brine_mixer)

    for i in m.fs.trains:
        build_ro_train(
            m.fs.train[i],
            prop_package=m.fs.properties,
            num_stages=num_stages,
            file=file,
        )

    m.fs.total_ro_pump_power = Expression(
        expr=sum(m.fs.train[i].total_pump_power for i in m.fs.trains)
    )

    for i, outlet in enumerate(outlet_list, start=1):

        sep_out = m.fs.ro_feed_separator.find_component(f"{outlet}")
        perm_mix_in = m.fs.ro_product_mixer.find_component(f"perm_inlet{i}")
        brine_mix_in = m.fs.ro_brine_mixer.find_component(f"brine_inlet{i}")
        a = Arc(
            source=sep_out,
            destination=m.fs.train[i].feed.inlet,
        )
        m.fs.add_component(f"sep_to_train{i}", a)
        a = Arc(
            source=m.fs.train[i].product.outlet,
            destination=perm_mix_in,
        )
        m.fs.add_component(f"train{i}_to_perm_mix", a)
        a = Arc(
            source=m.fs.train[i].disposal.outlet,
            destination=brine_mix_in,
        )
        m.fs.add_component(f"train{i}_to_brine_mix", a)

    if m.standalone:
        m.fs.recovery_vol_ro = Expression(
            expr=(m.fs.product.properties[0].flow_vol_phase["Liq"])
            / (m.fs.feed.properties[0].flow_vol_phase["Liq"])
        )

        m.fs.ro_feed_to_separator = Arc(
            source=m.fs.feed.outlet,
            destination=m.fs.ro_feed_separator.inlet,
        )

        m.fs.product_mixer_to_product = Arc(
            source=m.fs.ro_product_mixer.outlet,
            destination=m.fs.product.inlet,
        )

        m.fs.brine_mixer_to_disposal = Arc(
            source=m.fs.ro_brine_mixer.outlet,
            destination=m.fs.disposal.inlet,
        )

        TransformationFactory("network.expand_arcs").apply_to(m)

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )
    else:
        m.fs.recovery_vol_ro = Expression(
            expr=(m.fs.ro_product_mixer.mixed_state[0].flow_vol_phase["Liq"])
            / (m.fs.ro_feed_separator.mixed_state[0].flow_vol_phase["Liq"])
        )

    return m


def set_inlet_conditions(m, Qin=2637, Cin=0.5):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): m.num_trains
            * Qin
            * pyunits.gallons
            / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + 27,
        },
        hold_state=True,
    )


def set_ro_system_scaling(m):
    for i in m.fs.trains:
        set_ro_train_scaling(m.fs.train[i])


def set_ro_system_op_conditions(m):

    for i in m.fs.trains:
        set_ro_train_op_conditions(m.fs.train[i])
        if i != m.fs.trains.first():
            m.fs.ro_feed_separator.split_fraction[0, f"train{i}", "H2O"].fix(
                m.fs.ro_feed_separator.even_split
            )
            m.fs.ro_feed_separator.split_fraction[0, f"train{i}", "NaCl"].fix(
                m.fs.ro_feed_separator.even_split
            )
        else:
            m.fs.ro_feed_separator.split_fraction[0, f"train{i}", "H2O"].set_value(
                m.fs.ro_feed_separator.even_split
            )
            m.fs.ro_feed_separator.split_fraction[0, f"train{i}", "NaCl"].set_value(
                m.fs.ro_feed_separator.even_split
            )


def initialize_ro_system(m):

    if m.standalone:
        m.fs.feed.initialize()
        propagate_state(m.fs.ro_feed_to_separator)

    m.fs.ro_feed_separator.initialize()

    for i in m.fs.trains:
        a = m.fs.find_component(f"sep_to_train{i}")
        propagate_state(a)
        initialize_ro_train(m.fs.train[i])

        a = m.fs.find_component(f"train{i}_to_perm_mix")
        propagate_state(a)
        a = m.fs.find_component(f"train{i}_to_brine_mix")
        propagate_state(a)
        print(f"Initialized RO Train {i}")

    m.fs.ro_product_mixer.initialize()
    m.fs.ro_brine_mixer.initialize()

    if m.standalone:
        propagate_state(m.fs.product_mixer_to_product)
        m.fs.product.initialize()

        propagate_state(m.fs.brine_mixer_to_disposal)
        m.fs.disposal.initialize()


def add_ro_system_costing(m, costing_package=None, cost_RO=False):
    if costing_package is None:
        # costing_package = m.fs.costing
        m.fs.costing = costing_package = WaterTAPCosting()

    for i in m.fs.trains:
        add_ro_train_costing(m.fs.train[i], costing_package=costing_package, cost_RO=cost_RO)


def report_ro_system(m, w=30, add_costing=True):

    for i in m.fs.trains:
        # title = f"RO Train {i} Report"
        # side = int(((3 * w) - len(title)) / 2) - 1
        # header = "*" * side + f" {title} " + "*" * side
        # print(f"\n{header}\n")
        report_ro_train(m.fs.train[i], train_num=i, w=w, add_costing=add_costing)

    title = f"Overall RO System Performance"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "_" * side + f" {title} " + "_" * side
    print(f"\n\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    if m.standalone:
        print(
            f'{f"Total Perm Flow":<{w}s}{value(pyunits.convert(m.fs.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"Final Perm Conc":<{w}s}{value(pyunits.convert(m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
        print(
            f'{f"Total Brine Flow":<{w}s}{value(pyunits.convert(m.fs.disposal.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"Final Brine Conc":<{w}s}{value(pyunits.convert(m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
    else:
        print(
            f'{f"Total Perm Flow":<{w}s}{value(pyunits.convert(m.fs.ro_product_mixer.mixed_state[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"Final Perm Conc":<{w}s}{value(pyunits.convert(m.fs.ro_product_mixer.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
        print(
            f'{f"Total Brine Flow":<{w}s}{value(pyunits.convert(m.fs.ro_brine_mixer.mixed_state[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"Final Brine Conc":<{w}s}{value(pyunits.convert(m.fs.ro_brine_mixer.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
    print(f'{f"Overall Recovery":<{w}s}{value(m.fs.recovery_vol_ro)*100:<{w}.3f}{"%"}')
    print(
        f'{f"Total RO Pump Power":<{w}s}{value(pyunits.convert(m.fs.total_ro_pump_power, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )


def report_ro_system_pumps(m, w=30, add_costing=True):

    for i in m.fs.trains:
        for j in m.fs.train[i].stages:
            pump = m.fs.train[i].stage[j].pump
            title = f"Train {i} Stage {j} Pump Report"
            side = int(((3 * w) - len(title)) / 2) - 1
            header = "*" * side + f" {title} " + "*" * side
            print(f"\n{header}\n")
            report_pump(pump, w=w, add_costing=add_costing)


def main(add_costing=False):

    m = build_ro_system(num_trains=1, num_stages=2)
    set_ro_system_scaling(m)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=2637, Cin=0.5)
    set_ro_system_op_conditions(m)
    assert degrees_of_freedom(m) == 0
    initialize_ro_system(m)

    if add_costing:
        add_ro_system_costing(m)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_vol_phase["Liq"],
            name="SEC",
        )
        m.fs.costing.initialize()

    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_ro_system(m, w=40, add_costing=add_costing)
    # report_ro_system_pumps(m, w=20, add_costing=add_costing)
    return m


if __name__ == "__main__":
    m = main(add_costing=True)
