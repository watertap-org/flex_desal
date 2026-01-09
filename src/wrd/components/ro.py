from copy import deepcopy
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    TransformationFactory,
    value,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed, Product
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.models.unit_models import StateJunction
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)

from watertap.costing import WaterTAPCosting
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    PressureChangeType,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from wrd.utilities import load_config, get_config_value, get_config_file
from srp.utils import touch_flow_and_conc

default_ro_config = dict(
    has_pressure_change=True,
    pressure_change_type=PressureChangeType.fixed_per_stage,
    mass_transfer_coefficient=MassTransferCoefficient.calculated,
    concentration_polarization_type=ConcentrationPolarizationType.calculated,
    transformation_scheme="BACKWARD",
    transformation_method="dae.finite_difference",
    module_type="spiral_wound",
    finite_elements=10, #Change this to 1 to apply SD model equations to whole element
    has_full_reporting=True,
)


__all__ = [
    "build_ro",
    "initialize_ro",
    "set_ro_op_conditions",
    "set_ro_scaling",
    "report_ro",
    "add_ro_costing",
]

solver = get_solver()


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.ro = FlowsheetBlock(dynamic=False)

    build_ro(m.fs.ro, prop_package=m.fs.properties)

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)
    m.fs.brine = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.brine)

    m.fs.feed_to_ro = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.ro.feed.inlet,
    )
    m.fs.ro_to_product = Arc(
        source=m.fs.ro.product.outlet,
        destination=m.fs.product.inlet,
    )
    m.fs.ro_to_brine = Arc(
        source=m.fs.ro.disposal.outlet,
        destination=m.fs.brine.inlet,
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def build_ro(
    blk, prop_package=None, config=None, stage_num=1, file="wrd_inputs_8_19_21.yaml"
):

    if config is None:
        config = deepcopy(default_ro_config)
    else:
        for k, v in default_ro_config.items():
            if k not in config.keys():
                config[k] = v

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.stage_num = stage_num
    blk.config_data = load_config(get_config_file(file))
    config["property_package"] = prop_package

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.unit = ReverseOsmosis1D(**config)

    blk.product = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.product)
    blk.disposal = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.disposal)

    blk.feed_to_unit = Arc(
        source=blk.feed.outlet,
        destination=blk.unit.inlet,
    )
    blk.unit_to_product = Arc(
        source=blk.unit.permeate,
        destination=blk.product.inlet,
    )
    blk.unit_to_disposal = Arc(
        source=blk.unit.retentate,
        destination=blk.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_inlet_conditions(m, Qin=2637, Cin=0.5, file="wrd_inputs_8_19_21.yaml"):

    config_data = load_config(get_config_file(file))

    Pout = get_config_value(config_data, "pump_outlet_pressure", "pumps", f"pump_1")

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): Qin * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): Cin * pyunits.g / pyunits.L,
            ("pressure", None): Pout,
            ("temperature", None): 273.15 + 27,
        },
        hold_state=True,
    )


def set_ro_scaling(blk):
    """
    Add scaling to the units in the RO system
    """

    set_scaling_factor(blk.unit.feed_side.length, 1e-1)
    set_scaling_factor(blk.unit.feed_side.width, 1e-3)
    set_scaling_factor(blk.unit.area, 1e-5)
    set_scaling_factor(blk.unit.feed_side.area, 1e-5)
    set_scaling_factor(blk.unit.feed_side.spacer_porosity, 1e-1)
    # set_scaling_factor(blk.unit.feed_side.channel_height, 1e-5)

    for i, x in blk.unit.feed_side.mass_transfer_term.items():
        if i[3] == "NaCl":
            set_scaling_factor(x, 1e4)
        else:
            set_scaling_factor(x, 1)
    constraint_scaling_transform(blk.unit.feed_side.eq_dh, 1e-5)
    constraint_scaling_transform(blk.unit.eq_area, 1e-5)
    for i, c in blk.unit.feed_side.eq_K.items():
        set_scaling_factor(c, 1e4)


def set_ro_op_conditions(blk):

    # Set RO configuration for each stage
    print(f"Setting RO {blk.stage_num} operating conditions")
    blk.unit.A_comp.fix(
        get_config_value(
            blk.config_data, "A_comp", "reverse_osmosis_1d", f"stage_{blk.stage_num}"
        )
    )
    blk.unit.B_comp.fix(
        get_config_value(
            blk.config_data, "B_comp", "reverse_osmosis_1d", f"stage_{blk.stage_num}"
        )
    )

    blk.unit.feed_side.channel_height.fix(
        get_config_value(
            blk.config_data,
            "channel_height",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
    )
    blk.unit.feed_side.spacer_porosity.fix(
        get_config_value(
            blk.config_data,
            "spacer_porosity",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
    )

    blk.unit.feed_side.length.fix(
        get_config_value(
            blk.config_data,
            "number_of_elements_per_vessel",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
        * get_config_value(
            blk.config_data,
            "element_length",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
    )

    blk.unit.area.setub(1e6)
    blk.unit.width.setub(1e5)

    blk.unit.area.fix(
        get_config_value(
            blk.config_data,
            "element_membrane_area",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
        * get_config_value(
            blk.config_data,
            "number_of_vessels",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
        * get_config_value(
            blk.config_data,
            "number_of_elements_per_vessel",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
    )

    blk.unit.deltaP.fix(
        get_config_value(
            blk.config_data,
            "pressure_drop",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
    )

    blk.unit.recovery_vol_phase[0, "Liq"].set_value(  # Note this is unfixed!
        get_config_value(
            blk.config_data,
            "water_recovery_mass_phase",
            "reverse_osmosis_1d",
            f"stage_{blk.stage_num}",
        )
    )

    blk.unit.mixed_permeate[0].pressure.fix(20 * pyunits.psi)


def initialize_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_ro)

    initialize_ro(m.fs.ro)

    propagate_state(m.fs.ro_to_product)
    m.fs.product.initialize()

    propagate_state(m.fs.ro_to_brine)
    m.fs.brine.initialize()


def initialize_ro(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_unit)

    relax_bounds_for_low_salinity_waters(blk.unit)
    blk.unit.initialize()

    propagate_state(blk.unit_to_product)
    blk.product.initialize()

    propagate_state(blk.unit_to_disposal)
    blk.disposal.initialize()


def relax_bounds_for_low_salinity_waters(blk):
    blk.feed_side.cp_modulus.setub(5)
    for e in blk.feed_side.K:
        blk.feed_side.K[e].setub(0.01)
        blk.feed_side.K[e].setlb(1e-7)

    for e in blk.feed_side.cp_modulus:
        blk.feed_side.cp_modulus[e].setlb(1e-5)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.recovery_mass_phase_comp[e].setlb(1e-9)
            blk.recovery_mass_phase_comp[e].setub(1e-1)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.flux_mass_phase_comp[e].setlb(1e-9)
            blk.flux_mass_phase_comp[e].setub(1e-1)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "H2O":
            blk.recovery_mass_phase_comp[e].setlb(1e-4)
            blk.recovery_mass_phase_comp[e].setub(0.999)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "H2O":
            blk.flux_mass_phase_comp[e].setlb(1e-5)
            blk.flux_mass_phase_comp[e].setub(0.999)


def report_ro(blk, w=30):
    title = "RO Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    print(
        f'{f"Inlet Flow":<{w}s}{value(pyunits.convert(blk.feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallon / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Inlet Conc.":<{w}s}{value(pyunits.convert(blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.liter)):<{w}.3f}{"mg/L"}'
    )
    print(
        f'{f"Brine Flow":<{w}s}{value(pyunits.convert(blk.disposal.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallon / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Product Flow":<{w}s}{value(pyunits.convert(blk.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallon / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Recovery":<{w}s}{value(blk.unit.recovery_vol_phase[0, "Liq"])*100:<{w}.3f}{"%"}'
    )
    print(
        f'{f"Rejection":<{w}s}{value(blk.unit.rejection_phase_comp[0, "Liq", "NaCl"])*100:<{w}.3f}{"%"}'
    )
    print(
        f'{f"Perm Conc":<{w}s}{value(pyunits.convert(blk.unit.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.liter)):<{w}.3f}{f"mg/L"}'
    )
    print(
        f'{f"Brine Conc":<{w}s}{value(pyunits.convert(blk.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.liter)):<{w}.3f}{f"mg/L"}'
    )
    print(
        f'{f"∆P":<{w}s}{value(pyunits.convert(blk.unit.deltaP[0], to_units=pyunits.psi)):<{w}.3f}{f"psi"}'
    )
    print(
        f'{f"Membrane Area":<{w}s}{value(blk.unit.area):<{w}.3f}{f"{pyunits.get_units(blk.unit.area)}"}'
    )
    print(
        f'{f"Membrane Width":<{w}s}{value(blk.unit.width):<{w}.3f}{f"{pyunits.get_units(blk.unit.width)}"}'
    )
    print(
        f'{f"Membrane Length":<{w}s}{value(blk.unit.length):<{w}.3f}{f"{pyunits.get_units(blk.unit.length)}"}'
    )
    print(
        f'{f"Perm Backpressure":<{w}s}{value(pyunits.convert(blk.unit.mixed_permeate[0].pressure, to_units=pyunits.psi)):<{w}.3f}{f"psi"}'
    )


def add_ro_costing(blk, costing_package=None):

    if costing_package is None:
        m = blk.model()
        costing_package = m.fs.costing

    blk.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=costing_package)


def main(file="wrd_inputs_8_19_21.yaml"):

    m = build_system()
    set_ro_scaling(m.fs.ro)
    calculate_scaling_factors(m)
    set_inlet_conditions(m, file=file)
    set_ro_op_conditions(m.fs.ro)

    add_ro_costing(m.fs.ro)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])

    assert degrees_of_freedom(m) == 0
    initialize_system(m)

    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    report_ro(m.fs.ro)

    return m


if __name__ == "__main__":
    m = main()
