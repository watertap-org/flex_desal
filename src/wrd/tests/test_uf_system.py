import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.UF_system import main


@pytest.mark.component
def test_uf_system_8_19_21():
    m = main(num_trains=3, split_fraction=[0.385, 0.385, 0.23], Qin=10416, Cin=0.5)

    # Pump 1
    power = pyunits.convert(
        m.fs.uf_train[1].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 175 * pyunits.kW  # Measured value
    expected_power = 66 * pyunits.kW  # Modeled value
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kW

    # Pump 2
    power = pyunits.convert(
        m.fs.uf_train[2].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 175 * pyunits.kW  # Measured value
    expected_power = 66 * pyunits.kW  # Modeled value
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kW

    # Pump 3
    power = pyunits.convert(
        m.fs.uf_train[3].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 104 * pyunits.kW  # Measured value
    expected_power = 30 * pyunits.kW  # Modeled value
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kW

    total_power = pyunits.convert(m.fs.total_uf_pump_power, to_units=pyunits.kW)
    # expected_power = 454 * pyunits.kW # Measured value
    expected_power = 162 * pyunits.kW  # Modeled value
    assert pytest.approx(value(total_power), rel=0.15) == value(expected_power)  # kW


@pytest.mark.component
def test_uf_system_3_13_21():
    m = main(num_trains=3, split_fraction=[0.405, 0.405, 0.185], Qin=9764, Cin=0.5)

    # Pump 1
    power = pyunits.convert(
        m.fs.uf_train[1].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 178 * pyunits.kW  # Measured value
    expected_power = 66 * pyunits.kW  # Modeled value
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kW

    # Pump 2
    power = pyunits.convert(
        m.fs.uf_train[2].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 178 * pyunits.kW  # Measured value
    expected_power = 66 * pyunits.kW  # Modeled value
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kW

    # Pump 3
    power = pyunits.convert(
        m.fs.uf_train[3].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 79 * pyunits.kW  # Measured value
    expected_power = 26 * pyunits.kW  # Modeled value
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kW

    total_power = pyunits.convert(m.fs.total_uf_pump_power, to_units=pyunits.kW)
    # expected_power = 432 * pyunits.kW # Measured value
    expected_power = 154 * pyunits.kW  # Modeled value
    assert pytest.approx(value(total_power), rel=0.15) == value(expected_power)  # kW


@pytest.mark.component
def test_uf_system_with_costing():
    m = main(add_costing=True)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.129291  # kWh/m3


@pytest.mark.component
def test_uf_system_even_split():
    m = main()
    power = pyunits.convert(m.fs.total_uf_pump_power, to_units=pyunits.kW)
    expected_power = 76 * pyunits.kW  # Modeled Value
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kWh/m3
