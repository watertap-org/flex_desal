from pyomo.environ import (
    ConcreteModel,
    Expression,
    value,
    assert_optimal_termination,
    units as pyunits,
    value,
    Set,
    TransformationFactory,
)
from pyomo.network import Arc

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    Feed,
    Product,
    MomentumMixingType,
    Mixer,
    StateJunction,
)
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.costing import WaterTAPCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from wrd.utilities import load_config, get_config_file
from wrd.components.pump import *
from wrd.components.ro import *
from wrd.components.ro_stage import *
from srp.utils import touch_flow_and_conc

__all__ = [
    "build_ro_train",
    "set_ro_train_op_conditions",
    "set_ro_train_scaling",
    "initialize_ro_train",
    "report_ro_train",
    "add_ro_train_costing",
]

solver = get_solver()


def build_system(num_stages=3, file="wrd_ro_inputs_8_19_21.yaml"):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.ro_train = FlowsheetBlock(dynamic=False)
    build_ro_train(
        m.fs.ro_train, num_stages=num_stages, prop_package=m.fs.properties, file=file
    )

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)
    m.fs.brine = Product(property_package=m.fs.properties)

    # Arcs to connect the unit models
    m.fs.feed_to_train = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.ro_train.feed.inlet,
    )
    m.fs.train_to_product = Arc(
        source=m.fs.ro_train.product.outlet,
        destination=m.fs.product.inlet,
    )
    m.fs.train_to_brine = Arc(
        source=m.fs.ro_train.disposal.outlet,
        destination=m.fs.brine.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")  # changed from 1
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def set_inlet_conditions(m, Qin=2637, Cin=0.5, Tin=302, Pin=101325):

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): Pin,
            ("temperature", None): Tin,
        },
        hold_state=True,
    )


def build_ro_train(
    blk, num_stages=3, file="wrd_ro_inputs_8_19_21.yaml", prop_package=None
):

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.config_data = load_config(get_config_file(file))

    blk.stages = Set(initialize=range(1, num_stages + 1))
    blk.stage = FlowsheetBlock(blk.stages, dynamic=False)

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.mixer = Mixer(
        property_package=prop_package,
        inlet_list=[f"stage_{i}_to_product" for i in blk.stages],
        momentum_mixing_type=MomentumMixingType.none,
    )
    touch_flow_and_conc(blk.mixer)

    blk.product = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.product)
    blk.disposal = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.disposal)

    blk.recovery_vol = Expression(
        expr=blk.product.properties[0].flow_vol_phase["Liq"]
        / blk.feed.properties[0].flow_vol_phase["Liq"]
    )
    total_pump_power = 0

    for i in blk.stages:
        build_ro_stage(blk.stage[i], stage_num=i, file=file, prop_package=prop_package)
        total_pump_power += blk.stage[i].pump.unit.work_mechanical[0]

    for i in blk.stages:

        if i == blk.stages.first():
            ain = Arc(source=blk.feed.outlet, destination=blk.stage[i].feed.inlet)
            blk.add_component(f"feed_to_stage_{i}", ain)
            aout = Arc(
                source=blk.stage[i].disposal.outlet,
                destination=blk.stage[i + 1].feed.inlet,
            )
            blk.add_component(f"stage_{i}_to_stage_{i+1}", aout)

        elif i == blk.stages.last():
            aout_brine = Arc(
                source=blk.stage[i].disposal.outlet, destination=blk.disposal.inlet
            )
            blk.add_component(f"stage_{i}_to_brine", aout_brine)

        else:
            aout = Arc(
                source=blk.stage[i].disposal.outlet,
                destination=blk.stage[i + 1].feed.inlet,
            )
            blk.add_component(f"stage_{i}_to_stage_{i+1}", aout)

        mix_in = blk.mixer.find_component(f"stage_{i}_to_product")
        mix_arc = Arc(source=blk.stage[i].product.outlet, destination=mix_in)
        blk.add_component(f"stage_{i}_to_product", mix_arc)

    blk.mixer_to_product = Arc(source=blk.mixer.outlet, destination=blk.product.inlet)

    blk.total_pump_power = Expression(expr=total_pump_power)

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_ro_train_scaling(blk):

    for i in blk.stages:
        set_ro_stage_scaling(blk.stage[i])


def set_ro_train_op_conditions(blk):

    for i in blk.stages:
        set_ro_stage_op_conditions(blk.stage[i])

    blk.mixer.outlet.pressure[0].fix(101325)


def initialize_ro_train(blk):

    for i in blk.stages:
        if i == blk.stages.first():
            a = blk.find_component(f"feed_to_stage_{i}")
            propagate_state(a)
            initialize_ro_stage(blk.stage[i])
            a = blk.find_component(f"stage_{i}_to_stage_{i+1}")
            propagate_state(a)
        elif i == blk.stages.last():
            initialize_ro_stage(blk.stage[i])
            a = blk.find_component(f"stage_{i}_to_product")
            propagate_state(a)
            a = blk.find_component(f"stage_{i}_to_brine")
        else:
            initialize_ro_stage(blk.stage[i])
            a = blk.find_component(f"stage_{i}_to_stage_{i+1}")
            propagate_state(a)
        a = blk.find_component(f"stage_{i}_to_product")
        propagate_state(a)

    blk.mixer.initialize()
    propagate_state(blk.mixer_to_product)
    blk.product.initialize()


def initialize_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_train)

    initialize_ro_train(m.fs.ro_train)

    propagate_state(m.fs.train_to_product)
    m.fs.product.initialize()
    propagate_state(m.fs.train_to_brine)
    m.fs.brine.initialize()


def add_ro_train_costing(blk, costing_package=None):

    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    for i in blk.stages:
        add_ro_stage_costing(blk.stage[i], costing_package=costing_package)


def report_ro_train(blk, train_num=None, w=30):

    if train_num is None:
        title = "RO Train Report"
        title2 = f"Overall Train Performance"
    else:
        title = f"RO Train {train_num} Report"
        title2 = f"Overall Train {train_num} Performance"

    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    for i in blk.stages:
        title = f"Stage {i}"
        side = int(((3 * w) - len(title)) / 2) - 1
        header = "_" * side + f" {title} " + "_" * side
        print(f"\n\n{header}\n")
        print(
            f'{f"Stage {i} Feed Flow":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"Stage {i} Feed Conc.":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
        report_ro_stage(blk.stage[i], w=w)

    side = int(((3 * w) - len(title2)) / 2) - 1
    header = "(" * side + f" {title2} " + ")" * side
    print(f"\n\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{f"Total Pump Power":<{w}s}{value(pyunits.convert(blk.total_pump_power, to_units=pyunits.kilowatt)):<{w}.3f}{"kW"}'
    )

    print(f'{f"Overall Recovery":<{w}s}{value(blk.recovery_vol)*100:<{w}.3f}{"%"}')
    print(
        f'{f"Total Feed Flow":<{w}s}{value(pyunits.convert(blk.feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Inlet Feed Conc":<{w}s}{value(pyunits.convert(blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
    )
    print(
        f'{f"Total Perm Flow":<{w}s}{value(pyunits.convert(blk.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Final Perm Conc":<{w}s}{value(pyunits.convert(blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
    )
    print(
        f'{f"Total Brine Flow":<{w}s}{value(pyunits.convert(blk.disposal.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Final Brine Conc":<{w}s}{value(pyunits.convert(blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
    )
    print()
    for i, inlet in enumerate(blk.mixer.config.inlet_list, 1):
        sb = blk.mixer.find_component(f"{inlet}_state")
        print(
            f'{f"  Stage {i} Feed Flow":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"  Stage {i} Feed Conc":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )
        print(
            f'{f"  Stage {i} Perm Flow":<{w}s}{value(pyunits.convert(sb[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
        )
        print(
            f'{f"  Stage {i} Perm Conc":<{w}s}{value(pyunits.convert(sb[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
        )


def main(
    Qin=2637,
    Cin=0.528,
    Tin=302,
    Pin=101325,
    file="wrd_ro_inputs_8_19_21.yaml",
    add_costing=True,
):

    m = build_system(file=file)
    set_ro_train_scaling(m.fs.ro_train)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin, Tin=Tin, Pin=Pin)
    set_ro_train_op_conditions(m.fs.ro_train)

    initialize_system(m)
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    if add_costing:
        m.fs.costing = WaterTAPCosting()
        add_ro_train_costing(m.fs.ro_train)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_vol_phase["Liq"],
            name="SEC",
        )
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

    report_ro_train(m.fs.ro_train, w=30)

    return m


if __name__ == "__main__":
    m = main()

    # m = main(
    #     Qin=2452, Cin=0.503, Tin=295, Pin=101325, file="wrd_ro_inputs_3_13_21.yaml"
    # )
