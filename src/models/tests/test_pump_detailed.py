import pytest

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.control_volume_isothermal import ControlVolume0DBlock
from watertap.core.solvers import get_solver
from models.pump_detailed import Pump, VariableEfficiency

solver = get_solver()

# Build function with design flow and head as inputs
def build_pump_w_flow_head():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.unit = Pump(
        property_package=m.fs.properties,
        variable_efficiency=VariableEfficiency.Flow,
    )

    # Input flow and head
    feed_flow_vol = 0.126 * pyunits.m**3 / pyunits.s
    pump_head = 60.96 * pyunits.m
    density = 1000 * pyunits.kg / pyunits.m**3

    # Calculated feed conditions
    feed_flow_mass = feed_flow_vol * density
    feed_mass_frac_TDS = 0.035

    feed_pressure_in = 101325 * pyunits.Pa
    feed_pressure_out = feed_pressure_in + pump_head * density * 9.81 * pyunits.m / pyunits.s**2
    feed_temperature = 273.15 + 25

    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        feed_flow_mass * feed_mass_frac_TDS
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)

    m.fs.unit.system_curve_geometric_head.fix(4.57)
    m.fs.unit.ref_speed_fraction.fix(1.0)

    return m


# Build function with design flow and speed as inputs
def build_pump_w_flow_speed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.unit = Pump(
        property_package=m.fs.properties,
        variable_efficiency=VariableEfficiency.Flow,
    )

    # Input flow and speed
    feed_flow_vol = 0.126 * pyunits.m**3 / pyunits.s
    design_speed_fraction = 0.829
    density = 1000 * pyunits.kg / pyunits.m**3

    # Calculated feed conditions
    feed_flow_mass = feed_flow_vol * density
    feed_mass_frac_TDS = 0.035

    feed_pressure_in = 101325 * pyunits.Pa
    feed_temperature = 273.15 + 25

    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        feed_flow_mass * feed_mass_frac_TDS
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.design_speed_fraction.fix(design_speed_fraction)
    m.fs.unit.system_curve_geometric_head.fix(4.57)
    m.fs.unit.ref_speed_fraction.fix(1.0)

    return m


@pytest.mark.unit
def test_pump_w_flow_head():
    m = build_pump_w_flow_head()

    assert hasattr(m.fs.unit, "inlet")
    assert hasattr(m.fs.unit, "outlet")
    assert hasattr(m.fs.unit, "deltaP")  # this is just a reference
    assert hasattr(m.fs.unit.control_volume, "deltaP")

    m.fs.unit.initialize()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)

@pytest.mark.unit
def test_pump_w_flow_speed():
    m = build_pump_w_flow_speed()

    assert hasattr(m.fs.unit, "inlet")
    assert hasattr(m.fs.unit, "outlet")
    assert hasattr(m.fs.unit, "deltaP")  # this is just a reference
    assert hasattr(m.fs.unit.control_volume, "deltaP")

    m.fs.unit.initialize()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)

@pytest.mark.unit
def test_pump_w_head_speed():
    m = build_pump_w_flow_head()

    assert hasattr(m.fs.unit, "inlet")
    assert hasattr(m.fs.unit, "outlet")
    assert hasattr(m.fs.unit, "deltaP")  # this is just a reference
    assert hasattr(m.fs.unit.control_volume, "deltaP")

    m.fs.unit.initialize()
    assert degrees_of_freedom(m) == 0

    m.fs.unit.design_speed_fraction.fix(0.8184)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].unfix()
    m.fs.unit.control_volume.properties_in[0].mass_frac_phase_comp['Liq','TDS'].fix()

    results = solver.solve(m)
    assert_optimal_termination(results)