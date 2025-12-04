import pytest
from pyomo.environ import units as pyunits

from wrd.components.pump import main


@pytest.mark.component
def test_pump_PRO1_2_20():
    Qin = 2450 / 264.2 / 60  # gpm to m3/s
    Cin = 1082.1 * 0.5 / 1000  # us/cm to g/L
    # Pin = 101325 / 1e5  # Pa to bar # Concerned this is not correct! Other case gave 2x higher value and in psi
    Pin = 35.4 / 14.5  # psi to bar
    Pout = 153.5 / 14.5  # psi to bar
    expected_power = 201.7  # kW
    power = main(stage_num=1, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.1)


@pytest.mark.component
def test_pump_PRO1_8_19():
    # Conditions taken from ro_system_parameters_for_watertap.xlsx in Box August 19, 2021
    Qin = 2637 / 264.2 / 60  # gpm to m3/s
    Cin = 1055 * 0.5 / 1000  # us/cm to g/L
    Pin = 35.4 / 14.5  # psi to bar
    Pout = 141.9 / 14.5  # psi to bar
    expected_power = 196.25  # kW
    power = main(stage_num=1, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.1)


@pytest.mark.component
def test_pump_PRO1_3_13():
    Qin = 2452 / 264.2 / 60  # gpm to m3/s
    Cin = 1007 * 0.5 / 1000  # us/cm to g/L
    Pin = (
        35.46 / 14.5
    )  # Pa to bar # Concerned this is not correct! Other case gave 2x higher value and in psi
    Pout = 154 / 14.5  # psi to bar
    expected_power = 189.6  # kW
    power = main(stage_num=1, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.1)


@pytest.mark.component
def test_pump_PRO2_8_19():
    Qin = 1029 / 264.2 / 60  # gpm to m3/s
    Cin = 2496 * 0.5 / 1000  # us/cm to g/L
    Pin = (141.9 - 11.4) / 14.5  # psi to bar
    Pout = 160.5 / 14.5  # psi to bar
    expected_power = 22.7  # kW
    power = main(stage_num=2, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.1)


@pytest.mark.component
def test_pump_PRO2_2_20():
    Qin = 245.8 / 264.2 / 60  # gpm to m3/s
    Cin = 3055 * 0.5 / 1000  # us/cm to g/L
    Pin = (153.4 - 8.84) / 14.5  # psi to bar
    Pout = 157.8 / 14.5  # psi to bar
    expected_power = 9.7  # kW
    power = main(stage_num=2, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.15)


@pytest.mark.component
def test_pump_PRO2_3_13():
    Qin = 1041.4 / 264.2 / 60  # gpm to m3/s
    Cin = 2239 * 0.5 / 1000  # us/cm to g/L
    Pin = (154 - 9.9) / 14.5  # psi to bar
    Pout = 173.0 / 14.5  # psi to bar
    expected_power = 22.8  # kW
    power = main(stage_num=2, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.15)


# @pytest.mark.component
# def test_pump_TRO_2_20():
#     Qin = 123.3 / 264.2 / 60  # gpm to m3/s
#     Cin = 6346 * 0.5 / 1000  # us/cm to g/L
#     Pin = (157.8 - 6.54) / 14.5  # psi to bar
#     Pout = 156.8 / 14.5  # psi to bar
#     expected_power = 22.65  # kW
#     power = main(stage_num=3, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
#     assert power == pytest.approx(expected_power, rel=0.15)


# @pytest.mark.component
# def test_pump_TRO_8_19():
#     Qin = 383.6 / 264.2 / 60  # gpm to m3/s
#     Cin = 4847 * 0.5 / 1000  # us/cm to g/L
#     Pin = (160.5 - 7.2) / 14.5  # psi to bar
#     Pout = 164.4 / 14.5  # psi to bar
#     expected_power = 29.3  # kW
#     power = main(stage_num=3, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
#     assert power == pytest.approx(expected_power, rel=0.1)


# @pytest.mark.component
# def test_pump_TRO_3_13():
#     Qin = 517.3 / 264.2 / 60  # gpm to m3/s
#     Cin = 5540.58 * 0.5 / 1000  # us/cm to g/L
#     Pin = (172 - 7.0) / 14.5  # psi to bar
#     Pout = 173.0 / 14.5  # psi to bar
#     expected_power = 24.9  # kW
#     power = main(stage_num=3, Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
#     assert power == pytest.approx(expected_power, rel=0.15)
