import pytest
import pandas as pd
import os
from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
    value,
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
from models.pump_detailed import Pump, VariableEfficiency, PumpCurveDataType

solver = get_solver()

# Build function with design flow and head as inputs
def build_pump_w_flow_head():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.unit = Pump(
        property_package=m.fs.properties,
        variable_efficiency=VariableEfficiency.Flow,
        pump_curve_data_type=PumpCurveDataType.SurrogateCoefficent,
        head_surrogate_coeffs =  {0:114.22, 1:-410.6, 2:2729.2, 3:-8089.1},
        eff_surrogate_coeffs = {0:0.389, 1:-0.535, 2:41.373, 3:-138.82},
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
        pump_curve_data_type=PumpCurveDataType.SurrogateCoefficent,
        head_surrogate_coeffs =  {0:114.22, 1:-410.6, 2:2729.2, 3:-8089.1},
        eff_surrogate_coeffs = {0:0.389, 1:-0.535, 2:41.373, 3:-138.82},
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

# Three tests for different combinations of inputs
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

# Tests for take curves as dataset
@pytest.mark.unit
def test_data_points():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    pump_curves_filepath = os.path.join(os.path.dirname(__file__), "test_pump_curves_data.csv")

    m.fs.unit = Pump(
        property_package=m.fs.properties,
        variable_efficiency=VariableEfficiency.Flow,
        pump_curve_data_type=PumpCurveDataType.DataSet,
        # flow in m3/s and head in m
        pump_curves = pump_curves_filepath,
    )

    # Input flow and head
    feed_flow_vol = 0.1435 * pyunits.m**3 / pyunits.s
    pump_head = 87 * pyunits.m
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

    m.fs.unit.system_curve_geometric_head.fix(0)
    m.fs.unit.ref_speed_fraction.fix(1.0)

    assert hasattr(m.fs.unit, "inlet")
    assert hasattr(m.fs.unit, "surrogate_index")

    m.fs.unit.initialize()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)
    

# Test for different pumps
@pytest.mark.unit
def test_ro_feed_pump():
    # Test point from building out the pump curves?
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.unit = Pump(
        property_package=m.fs.properties,
        variable_efficiency=VariableEfficiency.Flow,
        pump_curve_data_type=PumpCurveDataType.SurrogateCoefficent,
        # pump_curves = os.path.join(os.path.dirname(__file__), "test_pump_curves_data.csv"),
        head_surrogate_coeffs =  {0:114.22, 1:-410.6, 2:2729.2, 3:-8089.1},
        eff_surrogate_coeffs = {0:0.389, 1:-0.535, 2:41.373, 3:-138.82},
    )
    
    # Input flow and head
    feed_flow_vol = 0.14 * pyunits.m**3 / pyunits.s
    pump_head = 172.3/3.28 * pyunits.m  
    density = 1000 * pyunits.kg / pyunits.m**3

    # Calculated feed conditions
    feed_flow_mass = feed_flow_vol * density
    feed_mass_frac_TDS = 0.035

    feed_pressure_in = 244074.4 * pyunits.Pa # 35.4 psi converted to Pa
    feed_pressure_out = feed_pressure_in + pump_head * density * (9.81 * pyunits.m / pyunits.s**2)
    feed_temperature = 273.15 + 25 # Why no units?

    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(feed_flow_mass * feed_mass_frac_TDS)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_flow_mass * feed_mass_frac_H2O)
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)

    m.fs.unit.system_curve_geometric_head.fix(0) # Could be a default value in model?
    m.fs.unit.ref_speed_fraction.fix(1.0) # Could be a default value in model?

    m.fs.unit.initialize()
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    assert value(m.fs.unit.ref_efficiency) == pytest.approx(0.825, abs=0.02)
    assert m.fs.unit.efficiency_pump[0].value == pytest.approx(0.76, abs=0.02) # 0.76 = 0.825*.95*.97


@pytest.mark.unit
def test_uf_pump():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.unit = Pump(
        property_package=m.fs.properties,
        variable_efficiency=VariableEfficiency.Flow,
        pump_curve_data_type=PumpCurveDataType.DataSet,
        pump_curves = os.path.join(os.path.dirname(__file__), "test_pump_curves_data_uf.csv"),
    )
    # Input flow and head
    feed_flow_vol = 0.142 * pyunits.m**3 / pyunits.s
    pump_head = 146.25/3.28 * pyunits.m
    density = 1000 * pyunits.kg / pyunits.m**3

    # Calculated feed conditions
    feed_flow_mass = feed_flow_vol * density
    feed_mass_frac_TDS = 0.035

    feed_pressure_in = 101325 * pyunits.Pa
    feed_pressure_out = feed_pressure_in + pump_head * density * (9.81 * pyunits.m / pyunits.s**2)
    feed_temperature = 273.15 + 25

    feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(feed_flow_mass * feed_mass_frac_TDS)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_flow_mass * feed_mass_frac_H2O)
    m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)

    m.fs.unit.system_curve_geometric_head.fix(4.57)
    m.fs.unit.ref_speed_fraction.fix(1.0)

    m.fs.unit.initialize()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)
    assert value(m.fs.unit.ref_efficiency) == pytest.approx(0.79, abs=0.02) #This doesn't account for geometric head
    assert m.fs.unit.efficiency_pump[0].value == pytest.approx(0.728, abs=0.02) # 0.728 = 0.79*.95*.97

if __name__ == "__main__":
    eff = test_data_points()
    print(eff)