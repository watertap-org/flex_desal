import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from wrd.components.pump import main


@pytest.mark.component
def test_pump_PRO_S1_8_19():
    # Conditions taken from ro_system_parameters_for_watertap.xlsx in Box August 19, 2021
    expected_power = 196.25 * pyunits.kW
    m = main(stage_num=1, date="8_19_21")
    power = pyunits.convert(
        m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
    )
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)


@pytest.mark.component
def test_pump_PRO_S2_8_19():
    expected_power = 22.7 * pyunits.kW
    expected_eta = 0.58
    m = main(stage_num=2, date="8_19_21")
    power = pyunits.convert(
        m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
    )
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.20)


@pytest.mark.component
def test_pump_TRO_8_19():
    expected_power = 29.3 * pyunits.kW # kW
    m = main(stage_num=3, date="8_19_21")
    power = pyunits.convert(
        m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
    )
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.20)


@pytest.mark.component
def test_pump_PRO_S1_3_13():
    expected_power = 189.6 * pyunits.kW 
    m = main(stage_num=1, date="3_13_21")
    power = pyunits.convert(
        m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
    )
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)


@pytest.mark.component
def test_pump_PRO_S2_3_13():
    extra_load = 3 * pyunits.kW
    expected_power = 22.8 * pyunits.kW - extra_load 
    m = main(stage_num=2, date="3_13_21")
    power = pyunits.convert(
        m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
    )
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.20)


@pytest.mark.component
def test_pump_TRO_3_13():
    expected_power = 24.9 * pyunits.kW 
    m = main(stage_num=3, date="2_20_21")
    power = pyunits.convert(
        m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
    )
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.20)


# This day is not being used for comparison for this report
# @pytest.mark.component
# def test_pump_PRO_S1_2_20():
#     expected_power = 201.7 * pyunits.kW
#     m = main(stage_num=3, date="2_20_21")
#     power = pyunits.convert(
#         m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
#     )
#     assert_units_consistent(power + expected_power)
#     assert value(power) == pytest.approx(value(expected_power), rel=0.15)


# @pytest.mark.component  # Operating at 80% speed instead of 100%
# def test_pump_PRO_S2_2_20():
#     expected_power = 9.7 * pyunits.kW
#     m = main(stage_num=3, date="2_20_21")
#     power = pyunits.convert(
#         m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
#     )
#     assert_units_consistent(power + expected_power)
#     assert value(power) == pytest.approx(value(expected_power), rel=0.20)


# @pytest.mark.component
# def test_pump_TRO_2_20():
#     expected_power = 22.65 * pyunits.kW
#     m = main(stage_num=3, date="2_20_21")
#     power = pyunits.convert(
#         m.fs.pump_system.pump.work_mechanical[0], to_units=pyunits.kW
#     )
#     assert_units_consistent(power + expected_power)
#     assert value(power) == pytest.approx(value(expected_power), rel=0.20)
