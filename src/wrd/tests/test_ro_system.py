import pytest

from wrd.components.ro_system import main


@pytest.mark.component
def test_ro_system():
    number_trains = 1
    Qin = 2637 / 264.2 / 60  # gpm to m3/s
    Cin = 1055 * 0.5 / 1000  # us/cm to g/L
    expected_power = 196.25 + 22.71 + 29.3
    expected_perm_flow_gpm = 1608 + 635 + 198
    power, perm_flow_gpm = main(number_trains, Qin, Cin)
    assert power == pytest.approx(expected_power, rel=0.5)
    assert perm_flow_gpm == pytest.approx(expected_perm_flow_gpm, rel=0.5)
