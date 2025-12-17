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

from wrd.components.uf_train import *
from wrd.components.pump import report_pump
from srp.utils import touch_flow_and_conc

solver = get_solver()


def build_uf_system(
    m=None,
    num_trains=3,
    split_fraction=None,
    prop_package=None,
    file="wrd_ro_inputs_8_19_21.yaml",
):

    if m is None:
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = NaClParameterBlock()
        m.fs.feed = Feed(property_package=m.fs.properties)
        touch_flow_and_conc(m.fs.feed)
        m.fs.product = Product(property_package=m.fs.properties)
        touch_flow_and_conc(m.fs.product)
        m.fs.disposal = Product(property_package=m.fs.properties)
        touch_flow_and_conc(m.fs.disposal)

    if prop_package is None:
        prop_package = m.fs.properties

    m.num_trains = num_trains
    m.fs.trains = Set(initialize=range(1, m.num_trains + 1))
    m.fs.train = FlowsheetBlock(m.fs.trains, dynamic=False)

    outlet_list = [f"train{i}" for i in m.fs.trains]

    m.fs.feed_separator = Separator(
        property_package=m.fs.properties,
        outlet_list=outlet_list,
        split_basis=SplittingType.componentFlow,
    )
    if split_fraction is None:
        # Even Split
        m.fs.feed_separator.feed_split = 1.0 / len(outlet_list)

    perm_inlet_list = [f"perm_inlet{i}" for i in m.fs.trains]

    m.fs.product_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=perm_inlet_list,
    )

    brine_inlet_list = [f"brine_inlet{i}" for i in m.fs.trains]

    m.fs.brine_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=brine_inlet_list,
    )

    # Could delete this tbh
    m.fs.recovery_vol = Expression(
        expr=(m.fs.product.properties[0].flow_vol_phase["Liq"])
        / (m.fs.feed.properties[0].flow_vol_phase["Liq"])
    )

    for i in m.fs.trains:
        build_uf_train(
            m.fs.train[i],
            prop_package=m.fs.properties,
            file=file,
        )

    for i, outlet in enumerate(outlet_list, start=1):

        sep_out = m.fs.feed_separator.find_component(f"{outlet}")
        perm_mix_in = m.fs.product_mixer.find_component(f"perm_inlet{i}")
        brine_mix_in = m.fs.brine_mixer.find_component(f"brine_inlet{i}")
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

    m.fs.feed_to_separator = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.feed_separator.inlet,
    )

    m.fs.product_mixer_to_product = Arc(
        source=m.fs.product_mixer.outlet,
        destination=m.fs.product.inlet,
    )

    m.fs.brine_mixer_to_disposal = Arc(
        source=m.fs.brine_mixer.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
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


def set_uf_system_scaling(m):

    for i in m.fs.trains:
        set_uf_train_scaling(m.fs.train[i])


def set_uf_system_op_conditions(m):

    for i in m.fs.trains:
        set_uf_train_op_conditions(m.fs.train[i])
        if i != m.fs.trains.first():
            m.fs.feed_separator.split_fraction[0, f"train{i}", "H2O"].fix(
                m.fs.feed_separator.feed_split
            )
            m.fs.feed_separator.split_fraction[0, f"train{i}", "NaCl"].fix(
                m.fs.feed_separator.feed_split
            )
        else:
            m.fs.feed_separator.split_fraction[0, f"train{i}", "H2O"].set_value(
                m.fs.feed_separator.feed_split
            )
            m.fs.feed_separator.split_fraction[0, f"train{i}", "NaCl"].set_value(
                m.fs.feed_separator.feed_split
            )

    m.fs.product_mixer.outlet.pressure[0].fix(101325)
    m.fs.brine_mixer.outlet.pressure[0].fix(101325)


def initialize_uf_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_separator)

    m.fs.feed_separator.initialize()

    for i in m.fs.trains:
        a = m.fs.find_component(f"sep_to_train{i}")
        propagate_state(a)
        initialize_uf_train(m.fs.train[i])

        a = m.fs.find_component(f"train{i}_to_perm_mix")
        propagate_state(a)
        a = m.fs.find_component(f"train{i}_to_brine_mix")
        propagate_state(a)

    m.fs.product_mixer.initialize()
    propagate_state(m.fs.product_mixer_to_product)
    m.fs.product.initialize()

    m.fs.brine_mixer.initialize()
    propagate_state(m.fs.brine_mixer_to_disposal)
    m.fs.disposal.initialize()


def add_uf_system_costing(m, costing_package=None):

    if costing_package is None:
        # costing_package = m.fs.costing
        m.fs.costing = costing_package = WaterTAPCosting()

    for i in m.fs.trains:
        add_uf_train_costing(m.fs.train[i], costing_package=costing_package)


def report_uf_system(m, w=30):

    for i in m.fs.trains:
        # title = f"RO Train {i} Report"
        # side = int(((3 * w) - len(title)) / 2) - 1
        # header = "*" * side + f" {title} " + "*" * side
        # print(f"\n{header}\n")
        report_uf_train(m.fs.train[i], train_num=i, w=w)

    title = f"Overall System Performance"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "_" * side + f" {title} " + "_" * side
    print(f"\n\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
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
    print(f'{f"Overall Recovery":<{w}s}{value(m.fs.recovery_vol)*100:<{w}.3f}{"%"}')


def report_uf_system_pumps(m, w=30):
    for i in m.fs.trains:
        pump = m.fs.train[i].pump
        title = f"Train {i} Pump Report"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "*" * side + f" {title} " + "*" * side
        print(f"\n{header}\n")
        report_pump(pump, w=w)


def main(add_costing=False):

    m = build_uf_system(num_trains=4)
    set_uf_system_scaling(m)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=2637, Cin=0.5)
    set_uf_system_op_conditions(m)
    assert degrees_of_freedom(m) == 0
    initialize_uf_system(m)

    if add_costing:
        add_uf_system_costing(m)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_vol_phase["Liq"],
            name="SEC",
        )
        m.fs.costing.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)
    report_uf_system(m, w=40)
    # report_uf_system_pumps(m, w=20)

    return m


if __name__ == "__main__":
    m = main(add_costing=True)
