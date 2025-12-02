from xml.parsers.expat import model
import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel,
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
from idaes.models.unit_models import Mixer, Separator, Product, Feed

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

import yaml
import os

solver = get_solver()


def build_inlet_stream(m):
    """Build the inlet stream for the RO system"""

    feed_conc = 2 * pyunits.g / pyunits.liter
    recovery = 0.15
    rho = 997.0 * pyunits.kg / pyunits.m**3

    perm_vol_flow = 9.1 * pyunits.m**3 / pyunits.day
    feed_vol_flow = pyunits.convert(
        perm_vol_flow / recovery, to_units=pyunits.m**3 / pyunits.s
    )
    perm_mass_flow = pyunits.convert(
        rho * perm_vol_flow, to_units=pyunits.kg / pyunits.s
    )

    feed_mass_flow_water = value(perm_mass_flow / recovery)
    feed_mass_flow_salt = value(
        pyunits.convert(feed_vol_flow * feed_conc, to_units=pyunits.kg / pyunits.s)
    )

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        feed_mass_flow_salt
    )  # 2 g/m^3 for NaCl
    m.fs.feed.properties[0].temperature.fix(298.15)  # 25 degrees Celsius
    m.fs.feed.properties[0].pressure.fix(101325)  # 1 atm


def set_operating_conditions(m):
    """Set the operating conditions for the RO system."""

    pressure = 225 * pyunits.psi
    # Pump 1 operating conditions
    m.fs.pump1.efficiency_pump.fix(0.8)
    m.fs.pump1.control_volume.properties_out[0].pressure.fix(pressure)

    mem_length = (1.016 - 2 * 0.0267) * pyunits.m
    mem_area = 7.2 * pyunits.m**2

    spacer_thickness = 34  # mil
    channel_height = spacer_thickness * 2.54e-5  # mil to m
    pressure_loss = -15 * pyunits.psi

    water_perm = 4.2e-12
    salt_perm = 3.5e-8
    porosity = 0.95

    m.fs.RO_stage_1.A_comp.fix(water_perm)
    m.fs.RO_stage_1.B_comp.fix(salt_perm)

    m.fs.RO_stage_1.feed_side.channel_height.fix(channel_height)
    m.fs.RO_stage_1.length.fix(mem_length)
    m.fs.RO_stage_1.feed_side.spacer_porosity.fix(porosity)
    # m.fs.RO_stage_1.area.fix(mem_area)

    # m.fs.RO_stage_1.feed_side.N_Re[0, 0].fix(400)
    # m.fs.RO_stage_1.deltaP.fix(pressure_loss)
    m.fs.RO_stage_1.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)  # 15% recovery
    # m.fs.RO_stage_1.permeate.pressure[0].fix(101325)  # 1 atm

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )


def add_connection(m):
    """Add connections between unit models in the RO system."""

    m.fs.feed_to_pump1 = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.pump1.inlet,
    )

    m.fs.pump1_to_RO_stage_1 = Arc(
        source=m.fs.pump1.outlet,
        destination=m.fs.RO_stage_1.inlet,
    )

    m.fs.RO_stage_1_to_permeate = Arc(
        source=m.fs.RO_stage_1.permeate,
        destination=m.fs.permeate.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def initialize_system(m):
    """Initialize the RO system."""

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_pump1)
    m.fs.pump1.initialize()

    propagate_state(m.fs.pump1_to_RO_stage_1)
    m.fs.RO_stage_1.initialize()

    propagate_state(m.fs.RO_stage_1_to_permeate)
    m.fs.permeate.initialize()


def build_system():
    """Build the system model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    # Feed stream to first pump
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.permeate = Product(property_package=m.fs.properties)

    # Feed pump to first stage RO
    m.fs.pump1 = Pump(property_package=m.fs.properties)

    # Inter stage pump between first and second stage RO
    # m.fs.pump2 = Pump(property_package=m.fs.properties)

    # Third stage pump to boost pressure before final RO
    # m.fs.pump3 = Pump(property_package=m.fs.properties)

    m.fs.RO_stage_1 = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        module_type="spiral_wound",
        finite_elements=10,
        has_full_reporting=True,
    )

    # m.fs.RO_stage_2 = ReverseOsmosis1D(
    #     property_package=m.fs.properties,
    #     has_pressure_change=True,
    #     pressure_change_type=PressureChangeType.calculated,
    #     mass_transfer_coefficient=MassTransferCoefficient.calculated,
    #     concentration_polarization_type=ConcentrationPolarizationType.calculated,
    #     transformation_scheme="BACKWARD",
    #     transformation_method="dae.finite_difference",
    #     module_type="spiral_wound",
    #     finite_elements=10,
    #     has_full_reporting=True,
    # )

    # m.fs.RO_stage_3 = ReverseOsmosis1D(
    #     property_package=m.fs.properties,
    #     has_pressure_change=True,
    #     pressure_change_type=PressureChangeType.calculated,
    #     mass_transfer_coefficient=MassTransferCoefficient.calculated,
    #     concentration_polarization_type=ConcentrationPolarizationType.calculated,
    #     transformation_scheme="BACKWARD",
    #     transformation_method="dae.finite_difference",
    #     module_type="spiral_wound",
    #     finite_elements=10,
    #     has_full_reporting=True,
    # )

    return m


def build_ro_stage(stage_number, m):
    """Build a single stage of the RO system."""

    # Set the stage number
    stage_name = f"stage_{stage_number}"

    # Create the RO stage
    ro_stage = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        module_type="spiral_wound",
        finite_elements=10,
        has_full_reporting=True,
    )

    # Set stage configuration
    # Get configuration from the config file
    (A_comp, B_comp) = (
        m.fs.config_data["reverse_osmosis_1d"]["stages"][stage_name]["A_comp"],
        m.fs.config_data["reverse_osmosis_1d"]["stages"][stage_name]["B_comp"],
    )

    return ro_stage


def build_wrd():
    """Build the WRD model with RO stages and pumps."""

    # Build system
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    # Feed stream to first pump and system permeate
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.permeate = Product(property_package=m.fs.properties)

    # Feed pump to first stage RO
    m.fs.pump1 = Pump(property_package=m.fs.properties)

    # First stage RO
    m.fs.RO_stage_1 = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        module_type="spiral_wound",
        finite_elements=10,
        has_full_reporting=True,
    )

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.7)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        0.0014
    )  # 2 g/m^3 for NaCl
    m.fs.feed.properties[0].temperature.fix(298.15)  # 25 degrees Celsius
    m.fs.feed.properties[0].pressure.fix(101325)

    pressure = 225 * pyunits.psi  # 50e5

    m.fs.pump1.efficiency_pump.fix(0.8)
    m.fs.pump1.control_volume.properties_out[0].pressure.fix(pressure)

    water_perm = 4.2e-12
    salt_perm = 3.5e-8
    porosity = 0.95
    number_of_elements = 7
    mem_length = (1.016 - 2 * 0.0267) * pyunits.m
    spacer_thickness = 34  # mil
    channel_height = spacer_thickness * 2.54e-5  # mil to m

    m.fs.RO_stage_1.A_comp.fix(water_perm)
    m.fs.RO_stage_1.B_comp.fix(salt_perm)

    m.fs.RO_stage_1.feed_side.channel_height.fix(channel_height)
    m.fs.RO_stage_1.length.fix(mem_length)
    m.fs.RO_stage_1.feed_side.spacer_porosity.fix(porosity)

    m.fs.RO_stage_1.recovery_vol_phase[0, "Liq"].fix(0.15)  # 15% recovery
    m.fs.RO_stage_1.permeate.pressure[0].fix(101325)  # 1 atm

    m.fs.feed.initialize()
    m.fs.feed_to_pump1 = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.pump1.inlet,
    )

    m.fs.pump1_to_RO_stage_1 = Arc(
        source=m.fs.pump1.outlet,
        destination=m.fs.RO_stage_1.inlet,
    )

    m.fs.RO_stage_1_to_permeate = Arc(
        source=m.fs.RO_stage_1.permeate,
        destination=m.fs.permeate.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    print("Degrees of freedom before initialization:", degrees_of_freedom(m))

    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_pump1)
    m.fs.pump1.initialize()

    propagate_state(m.fs.pump1_to_RO_stage_1)

    m.fs.RO_stage_1.initialize()

    propagate_state(m.fs.RO_stage_1_to_permeate)
    m.fs.permeate.initialize()

    iscale.calculate_scaling_factors(m.fs.RO_stage_1)

    print("Degrees of freedom before initialization:", degrees_of_freedom(m))

    print("Degrees of freedom:", degrees_of_freedom(m))
    results = solver.solve(m, tee=False)

    return m


def load_config(config):
    with open(config, "r") as file:
        return yaml.safe_load(file)


if __name__ == "__main__":

    m = build_wrd()
