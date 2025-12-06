import pytest
from wrd.components.pump import main


@pytest.mark.component
def test_pump_PRO_S1_8_19():
    # Conditions taken from ro_system_parameters_for_watertap.xlsx in Box August 19, 2021
    expected_power = 196.25  # kW
    power = main(stage_num=1, date="8_19_21")
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component
def test_pump_PRO_S2_8_19():
    expected_power = 22.7  # kW
    power = main(stage_num=2, date="8_19_21")
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component
def test_pump_TRO_8_19():
    expected_power = 29.3  # kW
    power = main(stage_num=3, date="8_19_21")
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component
def test_pump_PRO_S1_3_13():
    expected_power = 189.6  # kW
    power = main(stage_num=1, date="3_13_21")
    assert power == pytest.approx(expected_power, rel=0.15)


# @pytest.mark.component
def test_pump_PRO_S2_3_13():
    expected_power = 22.8  # kW
    power = main(stage_num=2, date="3_13_21")
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component
def test_pump_TRO_3_13():
    expected_power = 24.9  # kW
    power = main(stage_num=3, date="3_13_21")
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component
def test_pump_PRO_S1_2_20():
    expected_power = 201.7  # kW
    power = main(stage_num=1, date="2_20_21")
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component  # Operating at 80% speed instead of 100%
def test_pump_PRO_S2_2_20():
    expected_power = 9.7  # kW
    power = main(stage_num=2, date="2_20_21")
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component
def test_pump_TRO_2_20():
    expected_power = 22.65  # kW
    power = main(stage_num=3, date="2_20_21")
    assert power == pytest.approx(expected_power, rel=0.15)
    