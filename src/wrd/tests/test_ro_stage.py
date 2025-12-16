import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.ro_stage import main


@pytest.mark.component
def test_ro_stages_8_19_21():
    # Stage 1
    expected_power = 196.25 * pyunits.kW
    m = main(
        Qin=2637,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        stage_num=1,
        file="wrd_ro_inputs_8_19_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)

    # Stage 2
    # expected_power = 22.71 * pyunits.kW
    expected_power = 18 * pyunits.kW
    m = m = main(
        Qin=1029,
        Cin=1.2479,
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_ro_inputs_8_19_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)

    # Stage 3
    # expected_power = 29.3 * pyunits.kW
    expected_power = 34 * pyunits.kW
    m = main(
        Qin=384,
        Cin=4.847 / 2,
        Tin=302,
        Pin=(112.6 - 41.9) * pyunits.psi,
        stage_num=3,
        file="wrd_ro_inputs_8_19_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
