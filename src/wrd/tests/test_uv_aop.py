import pytest

from wrd.components.uv_aop import main


@pytest.mark.component
def test_uv_aop():
    m = main()
    assert (
        pytest.approx(m.fs.uv_aop_system.unit.power_consumption.value, rel=1e-4)
        == 98.329
    )
