from pyomo.environ import (
    ConcreteModel,
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
    StateJunction,
)
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.costing import WaterTAPCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from wrd.utilities import load_config, get_config_file
from wrd.components.pump import *
from wrd.components.UF_separator import *
from srp.utils import touch_flow_and_conc

__all__ = [
    "build_uf_train",
    "set_uf_train_op_conditions",
    "set_uf_train_scaling",
    "initialize_uf_train",
    "report_uf_train",
    "add_uf_train_costing",
]

solver = get_solver()


def build_system(file="wrd_uf_inputs_8_19_21.yaml"):
    # Will want to combine all inputs into one yaml instead of having separate ones
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.uf_train = FlowsheetBlock(dynamic=False)
    build_uf_train(m.fs.uf_train, prop_package=m.fs.properties, file=file)

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.uf_train = FlowsheetBlock(dynamic=False)
    build_uf_train(m.fs.uf_train, prop_package=m.fs.properties, file=file)

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)
    m.fs.disposal = Product(property_package=m.fs.properties)

    # Arcs to connect the unit models
    m.fs.feed_to_train = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.uf_train.feed.inlet,
    )
    m.fs.train_to_product = Arc(
        source=m.fs.uf_train.product.outlet,
        destination=m.fs.product.inlet,
    )
    m.fs.train_to_disposal = Arc(
        source=m.fs.uf_train.disposal.outlet,
        destination=m.fs.disposal.inlet,
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def build_uf_train(blk, file="wrd_uf_pump_inputs_8_19_21.yaml", prop_package=None):

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    name = (
        blk.name.split(".")[-1]
        .replace("_", " ")
        .replace("[", " ")
        .replace("]", "")
        .upper()
    )

    print(f'\n{f"=======> BUILDING {name} <=======":^60}\n')

    blk.config_data = load_config(get_config_file(file))

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.product = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.product)

    blk.disposal = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.disposal)

    blk.pump = FlowsheetBlock(dynamic=False)
    build_pump(blk.pump, file=file, prop_package=prop_package)
    blk.pump.config_data = (
        blk.config_data
    )  # Will need to revist how config data is being handled

    blk.UF = FlowsheetBlock(dynamic=False)
    build_separator(
        blk.UF, prop_package=prop_package, outlet_list=["product", "disposal"]
    )

    blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump.feed.inlet)
    blk.pump_to_UF = Arc(source=blk.pump.product.outlet, destination=blk.UF.feed.inlet)
    blk.UF_to_disposal = Arc(
        source=blk.UF.disposal.outlet, destination=blk.disposal.inlet
    )
    blk.UF_to_product = Arc(source=blk.UF.product.outlet, destination=blk.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_uf_train_scaling(blk):
    add_pump_scaling(blk.pump)
    # add_uf_scaling(blk.UF) # Seems like there are no variables to scale for separator


def set_inlet_conditions(m, Qin=2637, Cin=0.5, Tin=298, Pin=101325):
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): Pin,
            ("temperature", None): Tin,
        },
        hold_state=True,
    )


def set_uf_train_op_conditions(blk, split_fractions=None):
    set_pump_op_conditions(blk.pump)
    if split_fractions is None:
        split_fractions = {
            "product": {"H2O": 0.99, "NaCl": 0.99},
        }
    set_separator_op_conditions(blk.UF, split_fractions)


def initialize_system(m):
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_train)

    initialize_uf_train(m.fs.uf_train)

    propagate_state(m.fs.train_to_product)
    m.fs.product.initialize()
    propagate_state(m.fs.train_to_disposal)
    m.fs.disposal.initialize()


def initialize_uf_train(blk):
    blk.feed.initialize()

    propagate_state(blk.feed_to_pump)
    initialize_pump(blk.pump)

    propagate_state(blk.pump_to_UF)
    init_separator(blk.UF)

    propagate_state(blk.UF_to_disposal)
    blk.disposal.initialize()

    propagate_state(blk.UF_to_product)
    blk.product.initialize()


def add_uf_train_costing(blk, costing_package=None):

    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    add_pump_costing(blk.pump, costing_package=costing_package)
    # add_separator_costing(blk.ro,costing_package=costing_package) # Don't think there's anything to cost here required


def report_uf_train(blk, train_num=0, w=30):
    title = f"UF Train {train_num} Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(
        f'{f"Total Pump Power":<{w}s}{value(pyunits.convert(blk.pump.unit.work_mechanical[0], to_units=pyunits.kilowatt)):<{w}.3f}{"kW"}'
    )
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


def main(
    Qin=2637,
    Cin=0.528,
    Tin=302,
    Pin=101325,
    file="wrd_ro_inputs_8_19_21.yaml",
    add_costing=True,
):

    m = build_system(file=file)
    set_uf_train_scaling(m.fs.uf_train)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, Qin=Qin, Cin=Cin, Tin=Tin, Pin=Pin)
    set_uf_train_op_conditions(m.fs.uf_train)
    assert degrees_of_freedom(m) == 0
    initialize_system(m)
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    if add_costing:
        m.fs.costing = WaterTAPCosting()
        add_uf_train_costing(m.fs.uf_train)
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
    report_uf_train(m.fs.uf_train, w=30)

    return m


if __name__ == "__main__":
    m = main()
