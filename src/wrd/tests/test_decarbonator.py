import pytest

from wrd.components.decarbonator import main

@pytest.mark.component
def test_decarbonator():
    m = main()
    assert (
        pytest.approx(m.fs.decarb_system.unit.power_consumption.value, rel=1e-3)
        == 2.5
    )