import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from wrd.components.UF_train import main


@pytest.mark.component
def test_uf_train_with_costing():
    m = main(add_costing=True)


@pytest.mark.component
def test_uf_train_8_19_21():
    m = main(
        Qin=2637,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        file="wrd_ro_inputs_8_19_21.yaml",
        add_costing=True,
    )
    power = pyunits.convert(
        m.fs.uf_train.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 400 * pyunits.kW
    assert_units_consistent(power + expected_power)
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kWh/m3
