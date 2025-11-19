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
    list_badly_scaled_variables,
    extreme_jacobian_rows,
)

import yaml
import os


solver = get_solver()


def load_config(config):
    with open(config, "r") as file:
        return yaml.safe_load(file)


def get_config_value(   # Did not touch, don't think this needs to change - Jake
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
    #blk here is only the ro unit block
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
    m.fs.ro_system = FlowsheetBlock(dynamic=False)
    build_wrd_ro_system(m.fs.ro_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_wrd_ro_system(blk, prop_package=None, stage_num=1):
    """
    Build reverse osmosis system for WRD
    stage_num is the current stage number, which determines membrane properties 
    """
    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.ro_properties
    # Stage number
    blk.stage_num = stage_num

    # Feed stream, permeate, and brine 
    blk.feed = StateJunction(property_package=prop_package)
    blk.permeate = Product(property_package=prop_package)
    blk.brine = Product(property_package=prop_package)

    blk.recovery = Var( # Creating variable for RR. This may no longer be needed here, but moved to main flowsheet as an overall recovery from the different stages
        initialize=0.5,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Overall recovery",
    )

    @blk.Constraint(doc="Overall recovery constraint")
    def eq_recovery(b):
        return (
            b.recovery == b.permeate.properties[0].flow_vol_phase["Liq"] / b.feed.properties[0].flow_vol_phase["Liq"]
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
                # pressure_change_type=PressureChangeType.calculated, # Why this setting?
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
    blk.feed_to_RO = Arc(
        source=blk.feed.outlet, destination=blk.ro.inlet
    )

    blk.RO_to_permeate = Arc(
        source=blk.ro.outlet, destination=blk.permeate.inlet
    )
    blk.RO_to_brine = Arc(
        source=blk.ro.retentate, destination=blk.brine.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(blk)
    print("Degrees of freedom after adding units:", degrees_of_freedom(blk))


def set_inlet_conditions(blk, Qin=0.154, Cin=0.542):
    """
    Set the operation conditions for the RO system
    """
    Qin = (Qin) * pyunits.m**3 / pyunits.s  # Feed flow rate in m3/s
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    blk.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed.properties[0].pressure.fix(101325 * pyunits.Pa)  # 1 bar


def set_ro_system_op_conditions(blk):
    """
    Set the operation conditions for the RO system
    """
    s = blk.stage_num
    # Set RO configuration for each stage
    blk.ro.A_comp.fix(
        get_config_value(
            blk.config_data, "A_comp", "reverse_osmosis_1d", f"stage_{s}"
        )
    )
    blk.ro.B_comp.fix(
        get_config_value(
            blk.config_data, "B_comp", "reverse_osmosis_1d", f"stage_{s}"
        )
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
    blk.ro.e.deltaP.fix(
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
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    # RO Variables
    set_scaling_factor(blk.ro.feed_side.length, 1e-1)
    set_scaling_factor(blk.ro.feed_side.width, 1e-3)
    set_scaling_factor(blk.ro.area, 1e-5)
    set_scaling_factor(blk.ro.feed_side.spacer_porosity, 1e-1)
    # set_scaling_factor(blk.feed_side.channel_height, 1e-5)
    for i, x in blk.ro.feed_side.mass_transfer_term.items():
        if i[3] == "NaCl":
            set_scaling_factor(x, 1e4)
        else:
            set_scaling_factor(x, 1)
    constraint_scaling_transform(blk.feed_side.eq_dh, 1e-5)
    constraint_scaling_transform(blk.eq_area, 1e-5)
    for i, c in blk.ro.feed_side.eq_K.items():
        set_scaling_factor(c, 1e4)
    # constraint_scaling_transform(blk.feed_side.eq_K, 1e4)

    calculate_scaling_factors(blk)


def initialize_ro_system(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_RO)
    
    relax_bounds_for_low_salinity_waters(blk.ro)
    blk.ro.initialize()

    propagate_state(blk.RO_to_permeate)
    blk.permeate.initialize()

    propagate_state(blk.RO_to_brine)
    blk.brine.initialize()


def report_ro(blk, w=30): # This is not super informative yet
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
    

if __name__ == "__main__":
    m = build_system(number_trains=4, number_stages=3)
    set_inlet_conditions(m.fs.ro_system, Qin=4*0.154, Cin=0.542)
    set_ro_system_op_conditions(m.fs.ro_system)
    add_ro_scaling(m.fs.ro_system)
    initialize_ro_system(m.fs.ro_system)
    m.fs.obj = Objective(
        expr=m.fs.ro_system.permeate.properties[0].flow_vol_phase["Liq"]
    )
    results = solver.solve(m)
    assert_optimal_termination(results)

    print(f"{iscale.jacobian_cond(m.fs.ro_system):.2e}")
    # m.fs.ro_system.recovery.display()
    m.fs.ro_system.total_power_consumption.display()
    report_ro(m.fs.ro, w=40)