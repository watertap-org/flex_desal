import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from wrd.components.UF_train import main


@pytest.mark.component
def test_uf_train_with_costing():
    m = main(add_costing=True)


# Add tests for the full power UF pumps, half power UF pumps, and total flowrate divided by four.
@pytest.mark.component
def test_uf_train_3_13_21_full():
    m = main(
        Qin=3955,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        file="wrd_inputs_3_13_21.yaml",
        add_costing=True,
    )
    power = pyunits.convert(
        m.fs.uf_train.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 178 * pyunits.kW # Measured value
    expected_power = 63 * pyunits.kW  # Modeled value
    assert_units_consistent(power + expected_power)
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)


@pytest.mark.component
def test_uf_train_3_13_21_half():
    m = main(
        Qin=1785,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        file="wrd_inputs_3_13_21.yaml",
        add_costing=True,
    )
    power = pyunits.convert(
        m.fs.uf_train.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 79 * pyunits.kW # Measured value
    expected_power = 26 * pyunits.kW  # Modeled value
    assert_units_consistent(power + expected_power)
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)


@pytest.mark.component
def test_uf_train_3_13_21_one_train():
    # This test is using the average flowrate across four trains
    # i.e. 2452 gpm is total plant feed divided by 4
    # No actual UF pump sees this flowrate. The idea is to see the energy use of the uf system
    # can be apportioned to "one train"
    m = main(
        Qin=2452,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        file="wrd_inputs_3_13_21.yaml",
        add_costing=True,
    )
    power = pyunits.convert(
        m.fs.uf_train.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 108 * pyunits.kW # Total measured uf pump power divided by 4
    expected_power = 30 * pyunits.kW  # Modeled value
    assert_units_consistent(power + expected_power)
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)


# Just doing the divided train because technically these other flowrates are already tested in uf_system tests


@pytest.mark.component
def test_uf_train_8_19_21():
    m = main(
        Qin=2637,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        file="wrd_inputs_8_19_21.yaml",
        add_costing=True,
    )
    power = pyunits.convert(
        m.fs.uf_train.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    # expected_power = 114 * pyunits.kW  # Total measured uf pump power divided by 4
    expected_power = 32 * pyunits.kW  # Modeled value
    assert_units_consistent(power + expected_power)
    assert pytest.approx(value(power), rel=0.15) == value(expected_power)  # kWh/m3
