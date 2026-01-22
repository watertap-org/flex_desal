import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.ro_train import main


@pytest.mark.component
def test_ro_train_main_no_costing():
    m = main(add_costing=False)
    assert not hasattr(m.fs, "costing")


@pytest.mark.component
def test_ro_train1():
    # Input values from 3/13/21
    expected_power = 200.50 * pyunits.kW
    expected_product_flow = 2049.18 * pyunits.gal / pyunits.min
    expected_SEC = pyunits.convert(
        expected_power / expected_product_flow, to_units=pyunits.kWh / pyunits.m**3
    )

    m = main(
        Qin=2452,
        Cin=0.503,
        Tin=302,
        Pin=34.2 * pyunits.psi,
        file="wrd_inputs_3_13_21.yaml",
    )

    actual_power = pyunits.convert(m.fs.ro_train.total_pump_power, to_units=pyunits.kW)
    assert pytest.approx(value(actual_power), rel=1e-3) == value(expected_power)

    actual_product_flow = pyunits.convert(
        m.fs.ro_train.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )
    assert pytest.approx(value(actual_product_flow), rel=1e-3) == value(
        expected_product_flow
    )
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == value(expected_SEC)
