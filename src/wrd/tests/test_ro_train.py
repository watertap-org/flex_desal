import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.ro_train import main


@pytest.mark.component
def test_ro_train1_8_19_21():
    expected_power = (196.25 + 22.7 + 29.3) * pyunits.kW
    expected_product_flow = (1608 + 635 + 198) * pyunits.gal / pyunits.min
    expected_SEC = pyunits.convert(
        expected_power / expected_product_flow, to_units=pyunits.kWh / pyunits.m**3
    )

    m = main()

    actual_power = pyunits.convert(m.fs.ro_train.total_pump_power, to_units=pyunits.kW)
    assert pytest.approx(value(actual_power), rel=0.15) == value(expected_power)

    actual_product_flow = pyunits.convert(
        m.fs.ro_train.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )
    assert pytest.approx(value(actual_product_flow), rel=0.15) == value(
        expected_product_flow
    )
    assert pytest.approx(value(m.fs.costing.SEC), rel=0.15) == value(expected_SEC)


# NOTE: without the headloss prior to stage 3, this test will result
# in a negative deltaP for third stage pump and fail with costing
@pytest.mark.skip
def test_ro_train1_3_13_21():
    expected_power = (189.6 + 22.8 + 24.9) * pyunits.kW
    expected_product_flow = (1404.7 + 617 + 279) * pyunits.gal / pyunits.min
    expected_SEC = pyunits.convert(
        expected_power / expected_product_flow, to_units=pyunits.kWh / pyunits.m**3
    )

    m = main(Qin=2452, Cin=0.503, Tin=295, Pin=101325, file="wrd_inputs_3_13_21.yaml")

    actual_power = pyunits.convert(m.fs.ro_train.total_pump_power, to_units=pyunits.kW)
    assert pytest.approx(value(actual_power), rel=0.15) == value(expected_power)

    actual_product_flow = pyunits.convert(
        m.fs.ro_train.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )
    assert pytest.approx(value(actual_product_flow), rel=0.15) == value(
        expected_product_flow
    )
    assert pytest.approx(value(m.fs.costing.SEC), rel=0.15) == value(expected_SEC)


@pytest.mark.component
def test_ro_train_main_no_costing():
    m = main(add_costing=False)
    assert not hasattr(m.fs, "costing")
