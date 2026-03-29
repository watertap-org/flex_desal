import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.UF_system import main


@pytest.mark.component
def test_uf_system_with_costing():
    m = main(add_costing=True)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.1195  # kWh/m3


@pytest.mark.component
def test_uf_system_even_split():
    # Default tests 3 trains. No split fraction provided
    m = main()


@pytest.mark.component
def test_uf_system_three_trains():
    m = main(num_trains=3, split_fraction=[0.405, 0.405, 0.185], Qin=9764, Cin=0.528)

    # Pump 1
    power = pyunits.convert(
        m.fs.uf_train[1].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 107.99 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW

    # Pump 2
    power = pyunits.convert(
        m.fs.uf_train[2].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 106.4376 * pyunits.kW  # Idk why this is slightly diff from pump 1.
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW

    # Pump 3
    power = pyunits.convert(
        m.fs.uf_train[3].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 67.58 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW

    # Total
    total_power = pyunits.convert(m.fs.total_uf_pump_power, to_units=pyunits.kW)
    expected_power = 282 * pyunits.kW
    assert pytest.approx(value(total_power), rel=1e-3) == value(expected_power)  # kW


@pytest.mark.component
def test_uf_system_single_train():
    m = main(num_trains=1, split_fraction=[1.0], Qin=3955, Cin=0.528)
    # Pump 1
    power = pyunits.convert(
        m.fs.uf_train[1].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 106.455 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW


@pytest.mark.component
def test_uf_system_two_trains():
    m = main(num_trains=2, split_fraction=[0.6, 0.4], Qin=3955 / 0.6, Cin=0.528)
    # Pump 1
    power = pyunits.convert(
        m.fs.uf_train[1].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 106.455 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW

    # Pump 2
    power = pyunits.convert(
        m.fs.uf_train[2].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 78.28 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW


@pytest.mark.component
def test_uf_system_four_trains():
    m = main(num_trains=4, split_fraction=[0.3, 0.3, 0.3, 0.1], Qin=10000, Cin=0.528)
    # Pump 1
    power = pyunits.convert(
        m.fs.uf_train[1].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 84.13 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW

    # Pump 2
    power = pyunits.convert(
        m.fs.uf_train[2].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 84.13 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW

    # Pump 3
    power = pyunits.convert(
        m.fs.uf_train[3].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 84.13 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW

    # Pump 4
    power = pyunits.convert(
        m.fs.uf_train[4].pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    expected_power = 57.75 * pyunits.kW
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)  # kW
