import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from wrd.components.UF_train import main


@pytest.mark.component
def test_uf_train_with_costing():
    m = main(add_costing=True)


# Tests for full power UF pumps and half power UF pumps.
# Using 3/13/21 for input data at the moment
@pytest.mark.component
def test_uf_train_full():
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
    expected_power = 106.45 * pyunits.kW  # Modeled value
    assert_units_consistent(power + expected_power)
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)


@pytest.mark.component
def test_uf_train_half():
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

    expected_power = 67.33 * pyunits.kW  # Modeled value
    assert_units_consistent(power + expected_power)
    assert pytest.approx(value(power), rel=1e-3) == value(expected_power)


# Should decide whether uf pump testing should be in this file or in the test_pump file.
