import pytest
from pyomo.environ import units as pyunits

from wrd.components.pump import main


@pytest.mark.component
def test_pump():
    # Conditions taken from ro_system_parameters_for_watertap.xlsx in Box August 19, 2021
    Qin = 2637 / 264.2 / 60  # gpm to m3/s
    Cin = 1055 * 0.5 / 1000  # us/cm to g/L
    Pin = 35.4 / 14.5  # psi to bar
    Pout = 141.9 / 14.5  # psi to bar
    expected_power = 196.25  # kW
    power = main(Qin=Qin, Cin=Cin, Pin=Pin, Pout=Pout)
    assert power == pytest.approx(expected_power, rel=0.15)
