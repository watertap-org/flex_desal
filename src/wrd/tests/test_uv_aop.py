import pytest

from wrd.components.uv_aop import main
from pyomo.environ import value, units as pyunits


@pytest.mark.component
def test_uv_aop():
    m = main()
    expected_power = 98.32 * pyunits.kW  # Measured value
    assert pytest.approx(
        value(m.fs.uv_aop_system.unit.power_consumption.value), rel=1e-3
    ) == value(expected_power)
