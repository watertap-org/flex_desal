import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.UF_system import main


@pytest.mark.component
def test_uf_system_8_19_21():
    # Does this test actually pass any data specific to the date?
    m = main()
    power = pyunits.convert(m.fs.total_uf_pump_power, to_units=pyunits.kW)
    expected_power = 140 * pyunits.kW  # NOT VALUE FROM DATA
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kWh/m3




@pytest.mark.skip
def test_uf_system_with_costing():
    m = main()
    # Idk where the SEC values is coming from
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.33924  # kWh/m3
