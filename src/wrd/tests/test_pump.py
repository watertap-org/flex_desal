import pytest
from pyomo.environ import units as pyunits

from wrd.components.pump import main


@pytest.mark.component
def test_pump_PRO1_2_20():
    # Conditions taken from ro_system_parameters_for_watertap.xlsx in Box August 19, 2021
    Qin = 2450 / 264.2 / 60  # gpm to m3/s
    Cin = 1082.1 * 0.5 / 1000  # us/cm to g/L
    # Pin = 101325 / 1e5  # Pa to bar # Concerned this is not correct! Other case gave 2x higher value and in psi
    Pin = 35.4 / 14.5  # psi to bar
    Pout = 153.5 / 14.5  # psi to bar
    expected_power = 201.7  # kW
    power = main(Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.1)


@pytest.mark.component
def test_pump_PRO1_8_19():
    # Conditions taken from ro_system_parameters_for_watertap.xlsx in Box August 19, 2021
    Qin = 2637 / 264.2 / 60  # gpm to m3/s
    Cin = 1055 * 0.5 / 1000  # us/cm to g/L
    Pin = 35.4 / 14.5  # psi to bar
    Pout = 141.9 / 14.5  # psi to bar
    expected_power = 196.25  # kW
    power = main(Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.1)


@pytest.mark.component
def test_pump_PRO1_3_13():
    # Conditions taken from ro_system_parameters_for_watertap.xlsx in Box August 19, 2021
    Qin = 2452 / 264.2 / 60  # gpm to m3/s
    Cin = 1007 * 0.5 / 1000  # us/cm to g/L
    Pin = (
        35.46 / 14.5
    )  # Pa to bar # Concerned this is not correct! Other case gave 2x higher value and in psi
    Pout = 154 / 14.5  # psi to bar
    expected_power = 189.6  # kW
    power = main(Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.1)
