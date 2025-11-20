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

from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice

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


def get_config_value(
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


def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()

    m.fs.ro_train = FlowsheetBlock(dynamic=False)

    build_wrd_ro_system(m.fs.ro_train, prop_package=m.fs.ro_properties)

    return m


def build_wrd_ro_system(blk, prop_package=None):
    """
    Build reverse osmosis system for WRD
    """

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.ro_properties

    # Feed stream to first pump and system permeate
    blk.total_ro_feed = StateJunction(property_package=prop_package)
    blk.feed = StateJunction(property_package=prop_package)
    blk.permeate = Product(property_package=prop_package)
    blk.total_ro_product = StateJunction(property_package=prop_package)
    blk.no_of_stages = Param(
        initialize=3, mutable=True, doc="Number of RO stages in the system"
    )

    # WRD RO configurations input file. References to all values included in yml file
    # Get the absolute path of the current script
    current_script_path = os.path.abspath(__file__)
    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)
    # Get the parent directory of the current directory (one folder prior)
    parent_directory = os.path.dirname(current_directory)

    config = parent_directory + "/meta_data/wrd_ro_system_inputs.yaml"
    blk.config_data = load_config(config)

    add_ro_units(blk, prop_package=prop_package)
    print("Degrees of freedom after adding units:", degrees_of_freedom(blk))


def add_ro_units(blk, prop_package=None):

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.ro_properties

    # Feed pump to first stage RO
    blk.pump1 = Pump(property_package=prop_package)

    # Feed pump to second stage RO
    blk.pump2 = Pump(property_package=prop_package)

    # Feed pump to third stage RO
    blk.pump3 = Pump(property_package=prop_package)

    # Three stages of reverse osmosis
    for i in range(1, (blk.no_of_stages() + 1)):
        setattr(
            blk,
            f"ro_stage_{i}",
            ReverseOsmosis1D(
                property_package=prop_package,
                has_pressure_change=True,
                pressure_change_type=PressureChangeType.fixed_per_stage,
                mass_transfer_coefficient=MassTransferCoefficient.calculated,
                concentration_polarization_type=ConcentrationPolarizationType.calculated,
                transformation_scheme="BACKWARD",
                transformation_method="dae.finite_difference",
                module_type="spiral_wound",
                finite_elements=7,
                has_full_reporting=True,
            ),
        )

    # Add permeate mixer
    blk.permeate_mixer = Mixer(
        property_package=prop_package,
        inlet_list=[
            "ro_stage_1_permeate",
            "ro_stage_2_permeate",
            "ro_stage_3_permeate",
        ],
        energy_mixing_type=MixingType.extensive,
        momentum_mixing_type=MomentumMixingType.minimize,
    )

    add_ro_connections(blk)

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_ro_system_op_conditions(blk):
    """
    Set the operation conditions for the RO system
    """
    blk.total_ro_feed.inlet.temperature[0].fix(298.15)  # Fix feed temperature to 25 C
    blk.total_ro_feed.inlet.pressure[0].fix(4e5)  # Fix feed pressure to 4 bar

    feed_mass_flow_water = get_config_value(
        blk.config_data, "feed_flow_water", "feed_stream"
    ) * get_config_value(blk.config_data, "feed_density_water", "feed_stream")

    feed_mass_flow_salt = (
        get_config_value(blk.config_data, "feed_conductivity", "feed_stream")
        * get_config_value(
            blk.config_data, "feed_conductivity_conversion", "feed_stream"
        )
        * get_config_value(blk.config_data, "feed_flow_water", "feed_stream")
    )

    blk.total_ro_feed.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_mass_flow_water
    )
    blk.total_ro_feed.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_mass_flow_salt
    )


def build_ro_inlet_stream(blk, test=False):
    """Build the inlet stream for the RO system"""

    if test:
        # The feed stream is divided by the number of trains and vessels in stage 1
        equal_division_factor = get_config_value(
            blk.config_data, "number_of_trains", "reverse_osmosis_1d"
        )

        feed_mass_flow_water = (
            get_config_value(blk.config_data, "feed_flow_water", "feed_stream")
            * get_config_value(blk.config_data, "feed_density_water", "feed_stream")
            / equal_division_factor
        )

        feed_mass_flow_salt = (
            get_config_value(blk.config_data, "feed_conductivity", "feed_stream")
            * get_config_value(
                blk.config_data, "feed_conductivity_conversion", "feed_stream"
            )
            * get_config_value(blk.config_data, "feed_flow_water", "feed_stream")
            / equal_division_factor
        )

        feed_temperature = get_config_value(
            blk.config_data, "feed_temperature", "feed_stream"
        )
        feed_pressure = get_config_value(
            blk.config_data, "feed_pressure", "feed_stream"
        )

        blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            feed_mass_flow_water
        )
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            feed_mass_flow_salt
        )
        blk.feed.properties[0].temperature.fix(feed_temperature)
        blk.feed.properties[0].pressure.fix(feed_pressure)

    else:
        # Feed to RO is divided by number of trains
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            blk.total_ro_feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            / get_config_value(
                blk.config_data, "number_of_trains", "reverse_osmosis_1d"
            )
        )
        blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            blk.total_ro_feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
            / get_config_value(
                blk.config_data, "number_of_trains", "reverse_osmosis_1d"
            )
        )
        blk.feed.properties[0].temperature.fix(
            blk.total_ro_feed.properties[0].temperature
        )
        blk.feed.properties[0].pressure.fix(blk.total_ro_feed.properties[0].pressure)


def set_ro_operation_conditions(blk):
    """
    Set the operation conditions for the RO system
    """
    # Set pump operating conditions
    for i in range(1, (blk.no_of_stages() + 1)):
        pump = getattr(blk, f"pump{i}")

        pump.control_volume.properties_out[0].pressure.fix(
            get_config_value(
                blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{i}"
            )
        )

        pump.efficiency_pump.fix(
            get_config_value(blk.config_data, "pump_efficiency", "pumps", f"pump_{i}")
        )

    # Set RO configuration for each stage
    for i in range(1, (blk.no_of_stages() + 1)):
        ro_stage = getattr(blk, f"ro_stage_{i}")

        ro_stage.A_comp.fix(
            get_config_value(
                blk.config_data, "A_comp", "reverse_osmosis_1d", f"stage_{i}"
            )
        )
        ro_stage.B_comp.fix(
            get_config_value(
                blk.config_data, "B_comp", "reverse_osmosis_1d", f"stage_{i}"
            )
        )

        ro_stage.feed_side.channel_height.fix(
            get_config_value(
                blk.config_data, "channel_height", "reverse_osmosis_1d", f"stage_{i}"
            )
        )
        ro_stage.feed_side.spacer_porosity.fix(
            get_config_value(
                blk.config_data, "spacer_porosity", "reverse_osmosis_1d", f"stage_{i}"
            )
        )

        ro_stage.feed_side.length.fix(
            get_config_value(
                blk.config_data,
                "number_of_elements_per_vessel",
                "reverse_osmosis_1d",
                f"stage_{i}",
            )
            * get_config_value(
                blk.config_data, "element_length", "reverse_osmosis_1d", f"stage_{i}"
            )
        )

        ro_stage.area.setub(1e6)
        ro_stage.width.setub(1e5)

        ro_stage.area.fix(
            get_config_value(
                blk.config_data,
                "element_membrane_area",
                "reverse_osmosis_1d",
                f"stage_{i}",
            )
            * get_config_value(
                blk.config_data, "number_of_vessels", "reverse_osmosis_1d", f"stage_{i}"
            )
            * get_config_value(
                blk.config_data,
                "number_of_elements_per_vessel",
                "reverse_osmosis_1d",
                f"stage_{i}",
            )
        )

        ro_stage.deltaP.fix(
            get_config_value(
                blk.config_data, "pressure_drop", "reverse_osmosis_1d", f"stage_{i}"
            )
        )

        ro_stage.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(
            get_config_value(
                blk.config_data,
                "water_recovery_mass_phase",
                "reverse_osmosis_1d",
                f"stage_{i}",
            )
        )


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


def add_ro_connections(blk):
    """
    Add connections between the units in the RO system
    """

    # Connect feed to first pump
    blk.feed_to_pump1 = Arc(source=blk.feed.outlet, destination=blk.pump1.inlet)

    # Connect first pump to first RO stage
    blk.pump1_to_ro_stage_1 = Arc(
        source=blk.pump1.outlet, destination=blk.ro_stage_1.inlet
    )

    # Connect first RO stage to second pump
    blk.ro_stage_1_to_pump2 = Arc(
        source=blk.ro_stage_1.retentate, destination=blk.pump2.inlet
    )

    # Connect second pump to second RO stage
    blk.pump2_to_ro_stage_2 = Arc(
        source=blk.pump2.outlet, destination=blk.ro_stage_2.inlet
    )

    # Connect second RO stage to third pump
    blk.ro_stage_2_to_pump3 = Arc(
        source=blk.ro_stage_2.retentate, destination=blk.pump3.inlet
    )

    # Connect third pump to third RO stage
    blk.pump3_to_ro_stage_3 = Arc(
        source=blk.pump3.outlet, destination=blk.ro_stage_3.inlet
    )

    # Connect permeate from first and second stages to permeate mixer
    blk.ro_stage_1_to_permeate_mixer = Arc(
        source=blk.ro_stage_1.permeate,
        destination=blk.permeate_mixer.ro_stage_1_permeate,
    )

    blk.ro_stage_2_to_permeate_mixer = Arc(
        source=blk.ro_stage_2.permeate,
        destination=blk.permeate_mixer.ro_stage_2_permeate,
    )

    # Connect third RO stage to permeate mixer
    blk.ro_stage_3_to_permeate_mixer = Arc(
        source=blk.ro_stage_3.permeate,
        destination=blk.permeate_mixer.ro_stage_3_permeate,
    )

    # Connect permeate mixer to permeate product stream
    blk.permeate_mixer_to_permeate = Arc(
        source=blk.permeate_mixer.outlet, destination=blk.permeate.inlet
    )

    # Aggregate permeate from all trains to total product stream
    # Temperature constraint
    @blk.Constraint()
    def constrain_total_ro_product_temperature(blk):
        return (
            blk.total_ro_product.properties[0].temperature
            == blk.permeate.properties[0].temperature
        )

    @blk.Constraint()
    def constrain_total_ro_product_pressure(blk):
        return (
            blk.total_ro_product.properties[0].pressure
            == blk.permeate.properties[0].pressure
        )

    @blk.Constraint()
    def constrain_total_ro_product_flow_mass_phase_comp(blk):
        return blk.total_ro_product.properties[0].flow_vol_phase[
            "Liq"
        ] == blk.permeate.properties[0].flow_vol_phase["Liq"] * get_config_value(
            blk.config_data, "number_of_trains", "reverse_osmosis_1d"
        )

    # Fix total product TDS
    blk.total_ro_product.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0)


def add_ro_scaling(blk):
    """
    Add scaling to the units in the RO system
    """
    m = blk.model()

    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    for i in range(1, (blk.no_of_stages() + 1)):
        pump = getattr(blk, f"pump{i}")

        set_scaling_factor(pump.control_volume.work, 1e-3)

        ro_stage = getattr(blk, f"ro_stage_{i}")

        # Calculate RO scaling factors
        set_scaling_factor(ro_stage.feed_side.length, 1e-1)
        set_scaling_factor(ro_stage.feed_side.width, 1e-3)
        set_scaling_factor(ro_stage.area, 1e-5)
        set_scaling_factor(ro_stage.feed_side.spacer_porosity, 1e-1)
        set_scaling_factor(ro_stage.feed_side.channel_height, 1e-5)

        constraint_scaling_transform(ro_stage.feed_side.eq_dh, 1e-5)
        constraint_scaling_transform(ro_stage.eq_area, 1e-5)

        # for e in ro_stage.eq_flux_mass:
        #     if "NaCl" in e:
        #         print(f"Scaling NaCl flux constraint: {e}")
        #         constraint_scaling_transform(ro_stage.eq_flux_mass[e], 1e-33)
        #     else:
        #         # Different scaling for H2O or apply default scaling
        #         constraint_scaling_transform(ro_stage.eq_flux_mass[e], 1e-2)
        # for e in ro_stage.feed_side.eq_K:
        #     constraint_scaling_transform(ro_stage.feed_side.eq_K[e], 1e4)

    calculate_scaling_factors(blk)


def initialize_ro_units(blk):
    """
    Initialize the units in the RO system
    """
    calculate_scaling_factors(blk)

    # Initialize feed stream
    # blk.feed.initialize()
    print("Degrees of freedom after feed initialization:", degrees_of_freedom(blk))

    # Propagate feed state to pump
    propagate_state(blk.feed_to_pump1)

    # Initialize pumps and RO
    for i in range(1, (blk.no_of_stages() + 1)):
        pump = getattr(blk, f"pump{i}")
        pump.initialize()

        # Propagate state from pump to RO stage
        propagate_state(getattr(blk, f"pump{i}_to_ro_stage_{i}"))

        ro_stage = getattr(blk, f"ro_stage_{i}")

        print("Degrees of freedom after pump initialization:", degrees_of_freedom(blk))

        relax_bounds_for_low_salinity_waters(ro_stage)

        try:
            ro_stage.initialize()
        except:

            print_infeasible_constraints(ro_stage)
            # print(len(list_badly_scaled_variables(ro_stage)))
            # print("Degrees of freedom after RO initialization:", degrees_of_freedom(m))

            # Get badly scaled constraints (extreme row norms)
            # badly_scaled_constraints = extreme_jacobian_rows(
            #     m,  # your model
            #     scaled=True,  # use scaled Jacobian
            #     large=1e4,    # constraints with row norm >= this are considered large
            #     small=1e-4    # constraints with row norm <= this are considered small
            # )

            # # Print results
            # for norm, constraint in badly_scaled_constraints:
            #     print(f"Constraint {constraint.name}: row norm = {norm}")

        try:
            # Propagate state from RO stage to pump
            propagate_state(getattr(blk, f"ro_stage_{i}_to_pump{i+1}"))
        except:
            pass

        # Propagate state from RO stage to permeate mixer
        propagate_state(getattr(blk, f"ro_stage_{i}_to_permeate_mixer"))

    # Initialize permeate mixer
    blk.permeate_mixer.initialize()

    propagate_state(blk.permeate_mixer_to_permeate)
    blk.permeate.properties[0].flow_vol_phase
    blk.permeate.initialize()

    blk.total_ro_product.initialize()


if __name__ == "__main__":

    m = build_system()

    set_ro_system_op_conditions(m.fs.ro_train)
    m.fs.ro_train.total_ro_feed.initialize()

    build_ro_inlet_stream(m.fs.ro_train)
    print("Degrees of freedom after setting inlet stream:", degrees_of_freedom(m))

    set_ro_operation_conditions(m.fs.ro_train)
    print(
        "Degrees of freedom after setting operation conditions:", degrees_of_freedom(m)
    )

    add_ro_scaling(m.fs.ro_train)

    initialize_ro_units(m.fs.ro_train)

    results = solver.solve(m, tee=True)

    # m.fs.ro_train.total_ro_feed.display()
    # m.fs.ro_train.feed.display()

    print(f"{iscale.jacobian_cond(m.fs.ro_train):.2e}")

    assert False

    # try:
    #     results = solver.solve(m, tee=True)
    # except:
    #     print_infeasible_constraints(m)

    print("Degrees of freedom:", degrees_of_freedom(m))

    # Feed flowrate
    print("Feed \n", m.fs.ro_train.feed.properties[0].display())
    print()
    # Track flow rates and pressures for pumps and RO stages
    for i in range(1, (m.fs.ro_train.no_of_stages() + 1)):
        pump = getattr(m.fs.ro_train, f"pump{i}")
        print(f"Pump {i} - Inlet Pressure: {value(pump.inlet.pressure[0])} Pa")
        print(
            f"Pump {i} - Inlet Flow Rate: {value(pump.inlet.flow_mass_phase_comp[0,'Liq', 'H2O'])} kg/s"
        )
        print(f"Pump {i} - Outlet Pressure: {value(pump.outlet.pressure[0])} Pa")
        print(
            f"Pump {i} - Outlet Flow Rate: {value(pump.outlet.flow_mass_phase_comp[0,'Liq', 'H2O'])} kg/s"
        )
        ro_stage = getattr(m.fs.ro_train, f"ro_stage_{i}")
        print(f"RO Stage {i} - Inlet Pressure: {value(ro_stage.inlet.pressure[0])} Pa")
        print(
            f"RO Stage {i} - Inlet Flow Rate: {value(ro_stage.inlet.flow_mass_phase_comp[0,'Liq', 'H2O'])} kg/s"
        )
        print(
            f"RO Stage {i} - Retentate Flow: {value(ro_stage.retentate.flow_mass_phase_comp[0,'Liq', 'H2O'])} kg/s"
        )

    print("Permeate Stream \n", m.fs.ro_train.permeate.properties[0].display())

    print(
        "Permeate production from all trains:",
        m.fs.ro_train.permeate.properties[0].flow_vol_phase["Liq"]()
        * get_config_value(
            m.fs.ro_train.config_data, "number_of_trains", "reverse_osmosis_1d"
        ),
    )

    print(
        "Total Product Stream \n",
        m.fs.ro_train.total_ro_product.properties[0].display(),
    )
    print("Degrees of freedom:", degrees_of_freedom(m))
