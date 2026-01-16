from numpy import ones
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

from wrd.components.UF_train import *
from wrd.components.pump import report_pump
from srp.utils import touch_flow_and_conc

solver = get_solver()

__all__ = [
    "build_uf_system",
    "set_uf_system_scaling",
    "set_uf_system_op_conditions",
    "initialize_uf_system",
    "add_uf_system_costing",
    "report_uf_system",
    "report_uf_system_pumps",
]


def build_uf_system(
    m=None,
    num_trains=3,
    split_fraction=None,
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

    m.uf_num_trains = num_trains
    m.fs.uf_trains = Set(initialize=range(1, m.uf_num_trains + 1))
    m.fs.uf_train = FlowsheetBlock(m.fs.uf_trains, dynamic=False)

    outlet_list = [f"uf{i}" for i in m.fs.uf_trains]

    m.fs.uf_feed_separator = Separator(
        property_package=m.fs.properties,
        outlet_list=outlet_list,
        split_basis=SplittingType.componentFlow,
    )

    if split_fraction is None:
        # Even Split
        m.fs.uf_feed_separator.split_frac_input = (
            1.0 / len(outlet_list) * ones(len(outlet_list))
        )
    else:
        m.fs.uf_feed_separator.split_frac_input = split_fraction

    perm_inlet_list = [f"uf_prod_inlet{i}" for i in m.fs.uf_trains]

    m.fs.uf_product_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.minimize,
        inlet_list=perm_inlet_list,
    )

    brine_inlet_list = [f"uf_disp_inlet{i}" for i in m.fs.uf_trains]

    m.fs.uf_disposal_mixer = Mixer(
        property_package=m.fs.properties,
        momentum_mixing_type=MomentumMixingType.minimize,
        inlet_list=brine_inlet_list,
    )

    for i in m.fs.uf_trains:
        build_uf_train(
            m.fs.uf_train[i],
            prop_package=m.fs.properties,
            file=file,
        )

    m.fs.total_uf_pump_power = Expression(
        expr=sum(m.fs.uf_train[i].pump.unit.work_mechanical[0] for i in m.fs.uf_trains)
    )

    for i, outlet in enumerate(outlet_list, start=1):

        sep_out = m.fs.uf_feed_separator.find_component(f"{outlet}")
        perm_mix_in = m.fs.uf_product_mixer.find_component(f"uf_prod_inlet{i}")
        brine_mix_in = m.fs.uf_disposal_mixer.find_component(f"uf_disp_inlet{i}")
        a = Arc(
            source=sep_out,
            destination=m.fs.uf_train[i].feed.inlet,
        )
        m.fs.add_component(f"sep_to_uf{i}", a)
        a = Arc(
            source=m.fs.uf_train[i].product.outlet,
            destination=perm_mix_in,
        )
        m.fs.add_component(f"uf{i}_to_prod_mix", a)
        a = Arc(
            source=m.fs.uf_train[i].disposal.outlet,
            destination=brine_mix_in,
        )
        m.fs.add_component(f"uf{i}_to_disp_mix", a)

    if m.standalone:
        m.fs.uf_feed_to_separator = Arc(
            source=m.fs.feed.outlet,
            destination=m.fs.uf_feed_separator.inlet,
        )

        m.fs.recovery_vol_uf = Expression(
            expr=(m.fs.product.properties[0].flow_vol_phase["Liq"])
            / (m.fs.feed.properties[0].flow_vol_phase["Liq"])
        )
        m.fs.uf_product_mixer_to_product = Arc(
            source=m.fs.uf_product_mixer.outlet,
            destination=m.fs.product.inlet,
        )

        m.fs.uf_disposal_mixer_to_disposal = Arc(
            source=m.fs.uf_disposal_mixer.outlet,
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
        m.fs.recovery_vol_uf = Expression(
            expr=(m.fs.uf_product_mixer.mixed_state[0].flow_vol_phase["Liq"])
            / (m.fs.uf_feed_separator.mixed_state[0].flow_vol_phase["Liq"])
        )

    return m


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


def set_uf_system_scaling(m):

    for i in m.fs.uf_trains:
        set_uf_train_scaling(m.fs.uf_train[i])


def set_uf_system_op_conditions(m):

    for i in m.fs.uf_trains:
        set_uf_train_op_conditions(m.fs.uf_train[i])
        if i != m.fs.uf_trains.first():
            m.fs.uf_feed_separator.split_fraction[0, f"uf{i}", "H2O"].fix(
                m.fs.uf_feed_separator.split_frac_input[i - 1]
            )
            m.fs.uf_feed_separator.split_fraction[0, f"uf{i}", "NaCl"].fix(
                m.fs.uf_feed_separator.split_frac_input[i - 1]
            )
        else:
            m.fs.uf_feed_separator.split_fraction[0, f"uf{i}", "H2O"].set_value(
                m.fs.uf_feed_separator.split_frac_input[i - 1]
            )
            m.fs.uf_feed_separator.split_fraction[0, f"uf{i}", "NaCl"].set_value(
                m.fs.uf_feed_separator.split_frac_input[i - 1]
            )

    # m.fs.uf_product_mixer.outlet.pressure[0].fix(101325)
    # m.fs.uf_disposal_mixer.outlet.pressure[0].fix(101325)


def initialize_uf_system(m):

    if m.standalone:
        m.fs.feed.initialize()
        propagate_state(m.fs.uf_feed_to_separator)

    m.fs.uf_feed_separator.initialize()

    for i in m.fs.uf_trains:
        a = m.fs.find_component(f"sep_to_uf{i}")
        propagate_state(a)
        initialize_uf_train(m.fs.uf_train[i])

        a = m.fs.find_component(f"uf{i}_to_prod_mix")
        propagate_state(a)
        a = m.fs.find_component(f"uf{i}_to_disp_mix")
        propagate_state(a)

    m.fs.uf_product_mixer.initialize()
    m.fs.uf_disposal_mixer.initialize()

    if m.standalone:
        propagate_state(m.fs.uf_product_mixer_to_product)
        m.fs.product.initialize()

        propagate_state(m.fs.uf_disposal_mixer_to_disposal)
        m.fs.disposal.initialize()


def add_uf_system_costing(m, costing_package=None):

    if costing_package is None:
        # costing_package = m.fs.costing
        m.fs.costing = costing_package = WaterTAPCosting()

    for i in m.fs.uf_trains:
        add_uf_train_costing(m.fs.uf_train[i], costing_package=costing_package)


def report_uf_system(m, w=30):

    for i in m.fs.uf_trains:
        # title = f"RO Train {i} Report"
        # side = int(((3 * w) - len(title)) / 2) - 1
        # header = "*" * side + f" {title} " + "*" * side
        # print(f"\n{header}\n")
        report_uf_train(m.fs.uf_train[i], train_num=i, w=w)

    title = f"Overall UF System Performance"
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
            f'{f"Total Perm Flow":<{w}s}{value(pyunits.convert(m.fs.uf_product_mixer.mixed_state[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"Final Perm Conc":<{w}s}{value(pyunits.convert(m.fs.uf_product_mixer.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
        print(
            f'{f"Total Brine Flow":<{w}s}{value(pyunits.convert(m.fs.uf_disposal_mixer.mixed_state[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"Final Brine Conc":<{w}s}{value(pyunits.convert(m.fs.uf_disposal_mixer.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
    print(f'{f"Overall Recovery":<{w}s}{value(m.fs.recovery_vol_uf)*100:<{w}.3f}{"%"}')
    print(
        f'{f"Total UF Pump Power":<{w}s}{value(pyunits.convert(m.fs.total_uf_pump_power, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )


def report_uf_system_pumps(m, w=30):
    for i in m.fs.uf_trains:
        pump = m.fs.uf_train[i].pump
        title = f"UF Train {i} Pump Report"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "*" * side + f" {title} " + "*" * side
        print(f"\n{header}\n")
        report_pump(pump, w=w)


def main(
    add_costing=False,
    num_trains=3,
    split_fraction=None,
    Qin=10654,
    Cin=0.5,
    file="wrd_inputs_8_19_21.yaml",
):

    m = build_uf_system(num_trains=num_trains, split_fraction=split_fraction, file=file)
    set_uf_system_scaling(m)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin)
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
    report_uf_system(m)
    report_uf_system_pumps(m)
    return m


if __name__ == "__main__":
    m = main(
        num_trains=3,
        split_fraction=[0.4, 0.4, 0.2],
        Qin=10654,
        Cin=0.5,
        file="wrd_inputs_8_19_21.yaml",
    )
    # m = main(add_costing=False)
    m.fs.total_uf_pump_power.display()
    # m.fs.uf_feed_separator.display()
