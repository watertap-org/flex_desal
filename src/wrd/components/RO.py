from pyomo.environ import (
    ConcreteModel,
    Param,
    Objective,
    Expression,
    check_optimal_termination,
    value,
    assert_optimal_termination,
    units as pyunits,
    value,
    TransformationFactory,
)


from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import DiagnosticsToolbox
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    MixingType,
    MomentumMixingType,
    Mixer,
    Separator,
    Product,
    Feed,
    StateJunction,
)

from pyomo.network import Arc
import idaes.core.util.scaling as iscale

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    PressureChangeType,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
    get_scaling_factor,
    list_badly_scaled_variables,
    extreme_jacobian_rows,
)

import yaml
import os


solver = get_solver()


def load_config(config):
    with open(config, "r") as file:
        return yaml.safe_load(file)


def get_config_value(  # Did not touch, don't think this needs to change - Jake
    config,
    key,
    section,
    subsection=None,
):
    """
    Get a value from the configuration file.
    """

    if section in config:
        if subsection:
            if subsection in config[section]:
                if key in config[section][subsection]:
                    if (
                        isinstance(config[section][subsection][key], dict)
                        and "value" in config[section][subsection][key]
                        and "units" in config[section][subsection][key]
                    ):
                        return config[section][subsection][key]["value"] * getattr(
                            pyunits, config[section][subsection][key]["units"]
                        )
                    return config[section][subsection][key]
                else:
                    raise KeyError(
                        f"Key '{key}' not found in subsection '{subsection}' of section '{section}' of the configuration."
                    )
            else:
                raise KeyError(
                    f"Section '{section}' or subsection '{subsection}' not found in the configuration."
                )
        else:
            if key in config[section]:
                if (
                    isinstance(config[section][key], dict)
                    and "value" in config[section][key]
                    and "units" in config[section][key]
                ):
                    return config[section][key]["value"] * getattr(
                        pyunits, config[section][key]["units"]
                    )
                return config[section][key]
            else:
                raise KeyError(
                    f"Key '{key}' not found in section '{section}' of the configuration."
                )
    else:
        raise KeyError(f"Section '{section}' not found in the configuration.")


def relax_bounds_for_low_salinity_waters(blk):
    # blk here is only the ro unit block
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


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.ro_stage = FlowsheetBlock(dynamic=False)
    build_wrd_ro(m.fs.ro_stage, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_wrd_ro(
    blk, prop_package=None, stage_num=1
):  # blk should be flowsheet with pumps and ro stages
    """
    Build reverse osmosis system for WRD
    stage_num is the current stage number, which determines membrane properties
    """
    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.ro_properties
    # Stage number
    blk.stage_num = stage_num

    # Feed stream, permeate, and retentate
    blk.feed = StateJunction(property_package=prop_package)
    blk.permeate = StateJunction(property_package=prop_package)
    blk.retentate = StateJunction(property_package=prop_package)

    blk.recovery = Var(  # Creating variable for RR. This may no longer be needed here, but moved to main flowsheet as an overall recovery from the different stages
        initialize=0.5,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Overall recovery",
    )

    @blk.Constraint(doc="Overall recovery constraint")
    def eq_recovery(b):
        return (
            b.recovery
            == b.permeate.properties[0].flow_vol_phase["Liq"]
            / b.feed.properties[0].flow_vol_phase["Liq"]
        )

    # WRD RO configurations input file. References to all values included in yml file
    # Get the absolute path of the current script
    current_script_path = os.path.abspath(__file__)
    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)
    # Get the parent directory of the current directory (one folder prior)
    parent_directory = os.path.dirname(current_directory)

    config = parent_directory + "/meta_data/wrd_ro_inputs.yaml"
    blk.config_data = load_config(config)
    """
    add_ro_units(train, prop_package=prop_package)
    Add RO units to a single RO train
    """
    blk.ro = ReverseOsmosis1D(
        property_package=prop_package,
        has_pressure_change=True,
        #pressure_change_type=PressureChangeType.calculated, # Why this setting?
        pressure_change_type=PressureChangeType.fixed_per_stage,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        module_type="spiral_wound",
        finite_elements=7,
        has_full_reporting=True,
    )

    """
    add_ro_connections(train)
    Add connections between the units in the RO system
    """
    # Connect permeate mixer to permeate product stream
    blk.feed_to_RO = Arc(source=blk.feed.outlet, destination=blk.ro.inlet)
    blk.RO_to_permeate = Arc(source=blk.ro.permeate, destination=blk.permeate.inlet)
    blk.RO_to_retentate = Arc(source=blk.ro.retentate, destination=blk.retentate.inlet)
    TransformationFactory("network.expand_arcs").apply_to(blk)
    print("Degrees of freedom after adding units:", degrees_of_freedom(blk))


def set_inlet_conditions(blk, Qin=0.154, Cin=0.542, P_in=10.6):
    """
    Set the operation conditions for the RO stage
    """
    Qin = (Qin) * pyunits.m**3 / pyunits.s  # Feed flow rate in m3/s
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    blk.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed.properties[0].pressure.fix(P_in * pyunits.bar)
    blk.feed.properties[0].flow_vol  # Touching
    # m = blk.model()
    # m.fs.ro_properties.flow_mass_phase_comp["Liq","H2O"]
    # m.fs.ro_properties.flow_mass_phase_comp["Liq","NaCl"]

def set_ro_op_conditions(blk):
    """
    Set the operation conditions for the RO system
    """
    s = blk.stage_num
    # Set RO configuration for each stage
    blk.ro.A_comp.fix(
        get_config_value(blk.config_data, "A_comp", "reverse_osmosis_1d", f"stage_{s}")
    )
    blk.ro.B_comp.fix(
        get_config_value(blk.config_data, "B_comp", "reverse_osmosis_1d", f"stage_{s}")
    )
    blk.ro.feed_side.channel_height.fix(
        get_config_value(
            blk.config_data,
            "channel_height",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
    )
    blk.ro.feed_side.spacer_porosity.fix(
        get_config_value(
            blk.config_data,
            "spacer_porosity",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
    )
    blk.ro.feed_side.length.fix(
        get_config_value(
            blk.config_data,
            "number_of_elements_per_vessel",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
        * get_config_value(
            blk.config_data,
            "element_length",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
    )
    blk.ro.area.fix(
        get_config_value(
            blk.config_data,
            "element_membrane_area",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
        * get_config_value(
            blk.config_data,
            "number_of_vessels",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
        * get_config_value(
            blk.config_data,
            "number_of_elements_per_vessel",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
    )
    blk.ro.deltaP.fix(
        get_config_value(
            blk.config_data, "pressure_drop", "reverse_osmosis_1d", f"stage_{s}"
        )
    )
    blk.ro.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(
        get_config_value(
            blk.config_data,
            "water_recovery_mass_phase",
            "reverse_osmosis_1d",
            f"stage_{s}",
        )
    )


def add_ro_scaling(blk):
    """
    Add scaling to the units in the RO system
    """
    # Properties
    m = blk.model()
    # print(get_scaling_factor(m.fs.ro_stage.feed.properties[0].flow_mass_phase_comp["Liq","H2O"]))
    # m.fs.ro_stage.feed.properties[0].set_default_scaling(
    #     "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    # )
    # print(get_scaling_factor(m.fs.ro_stage.feed.properties[0].flow_mass_phase_comp["Liq","H2O"]))
    # m.fs.ro_stage.feed.properties[0].set_default_scaling(
    #     "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    # )
    # print(get_scaling_factor(m.fs.ro_stage.feed.properties[0].flow_mass_phase_comp["Liq","H2O"])) # Will not work until feed.initialize!
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    # print(get_scaling_factor(m.fs.ro_stage.feed.properties[0].flow_mass_phase_comp["Liq","H2O"]))
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    # RO Variables
    set_scaling_factor(blk.ro.feed_side.length, 1e-1)
    set_scaling_factor(blk.ro.feed_side.width, 1e-3)
    set_scaling_factor(blk.ro.area, 1e-5)
    set_scaling_factor(blk.ro.feed_side.area,1e-5) # Why is there area and feed_side.area?
    set_scaling_factor(blk.ro.feed_side.spacer_porosity, 1e-1)
    # set_scaling_factor(blk.feed_side.channel_height, 1e-5)
    for i, x in blk.ro.feed_side.mass_transfer_term.items():
        if i[3] == "NaCl":  
            set_scaling_factor(x, 1e4)
        else:
            set_scaling_factor(x, 1)
    constraint_scaling_transform(blk.ro.feed_side.eq_dh, 1e-5)
    constraint_scaling_transform(blk.ro.eq_area, 1e-5)
    for i, c in blk.ro.feed_side.eq_K.items():
        set_scaling_factor(c, 1e4)
    # constraint_scaling_transform(blk.feed_side.eq_K, 1e4)

    calculate_scaling_factors(m)


def initialize_ro(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_RO)

    relax_bounds_for_low_salinity_waters(blk.ro)
    dt = DiagnosticsToolbox(blk.ro)
    dt.report_structural_issues()
    print(blk.ro.width.bounds)
    print(blk.ro.area.value / blk.ro.length.value)
    # Is A / L > W_max ???
    # relax that constraint?
    blk.ro.width.bounds = (0.1,3000)
    blk.ro.initialize()

    propagate_state(blk.RO_to_permeate)
    blk.permeate.initialize()

    propagate_state(blk.RO_to_retentate)
    blk.retentate.initialize()


def report_ro(blk, w=30):  # This is not super informative yet
    title = "RO Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    total_flow = blk.feed.properties[0].flow_vol
    print(
        f'{f"Total Flow Rate (MGD)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.Mgallons /pyunits.day)):<{w}.3f}{"MGD"}'
    )
    print(f'{f"Total Flow Rate (m3/s)":<{w}s}{value(total_flow):<{w}.3e}{"m3/s"}')
    print(
        f'{f"Total Flow Rate (gpm)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )

def check_scaling(blk, descend_into=True, w=(60, 20, 20, 20)):
    def format_entries(name, unscaled_value, scaling_factor, scaled_value):

        if "e" in repr(unscaled_value) and unscaled_value is not None:
            unscaled_value = f"{unscaled_value:.4e}"
        else:
            if unscaled_value == 0:
                unscaled_value = "<zero>"
            elif abs(unscaled_value) < 5e-3:
                unscaled_value = f"{unscaled_value:.4e}"
            else:
                unscaled_value = f"{unscaled_value:.4f}"
        if scaling_factor is None:
            scaling_factor = "NONE"
        elif "e" in repr(scaling_factor) and scaling_factor is not None:
            scaling_factor = f"{scaling_factor:.4e}"
        else:
            if scaling_factor == 0:
                scaling_factor = "<zero>"
            elif abs(scaling_factor) < 5e-3:
                scaling_factor = f"{scaling_factor:.4e}"
            else:
                scaling_factor = f"{scaling_factor:.4f}"
        if "e" in repr(scaled_value) and scaled_value is not None:
            scaled_value = f"{scaled_value:.4e}"
        else:
            if scaled_value == 0:
                scaled_value = "<zero>"
            elif abs(scaled_value) < 5e-3:
                scaled_value = f"{scaled_value:.4e}"
            else:
                scaled_value = f"{scaled_value:.4f}"
        print(
            f"{name:<{w[0]}} | {unscaled_value:>{w[1]}} | {scaling_factor:>{w[2]}} | {scaled_value:>{w[3]}} |"
        )

    print(f"\n\n")
    header = f"{'VAR NAME':<{w[0]}} | {'UNSCALED VALUE':>{w[1]}} | {'SCALING FACTOR':>{w[2]}} | {'SCALED VALUE':>{w[3]}} |"
    print(header)
    print("=" * len(header))
    for v in blk.component_objects(Var, descend_into=descend_into):
        if "._" in v.name:
            continue
        if v.is_indexed():
            for i, vv in v.items():
                x = v.name.replace("[0.0]", "")
                x = x.split(".")
                vname = f"{x[-2]}.{x[-1]}[{i}]"
                sf = get_scaling_factor(vv)
                try:
                    var_val = value(vv)
                except ValueError:
                    var_val = None
                if var_val is None:
                    continue
                if sf is None:
                    sval = var_val
                    format_entries(vname, var_val, sf, sval)
                    continue
                sval = var_val * sf
                format_entries(vname, var_val, sf, sval)
        else:
            vname = v.name.split(".")[-1]

            sf = get_scaling_factor(v)
            try:
                var_val = value(v)
            except ValueError:
                var_val = None
            if var_val is None:
                continue
            if sf is None:
                sval = var_val
                format_entries(vname, var_val, sf, sval)
                continue
            sval = var_val * sf
            format_entries(vname, var_val, sf, sval)


if __name__ == "__main__":
    m = build_system()  # optional input of stage_num
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_inlet_conditions(
        m.fs.ro_stage, Qin=0.154, Cin=0.542, P_in=10.6
    )  # ro_system just adds feed,retentate,perm. May consider renaming to avoid implying pumps are included?
    set_ro_op_conditions(m.fs.ro_stage)
    print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    add_ro_scaling(m.fs.ro_stage)
    check_scaling(m.fs.ro_stage.ro)
    dt = DiagnosticsToolbox(m.fs.ro_stage)
    dt.report_structural_issues()
    initialize_ro(m.fs.ro_stage)
    m.fs.obj = Objective(
        expr=m.fs.ro_stage.permeate.properties[0].flow_vol_phase["Liq"]
    )
    results = solver.solve(m)
    assert_optimal_termination(results)

    print(f"{iscale.jacobian_cond(m.fs.ro_stage):.2e}")
    m.fs.ro_stage.recovery.display()
    report_ro(m.fs.ro_stage, w=40)
