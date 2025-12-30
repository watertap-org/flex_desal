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
        file="wrd_inputs_8_19_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)

    # Stage 2
    expected_power = 22.71 * pyunits.kW
    m = m = main(
        Qin=1175,  # 1029 -> measure value. 1175 is flow in modeled in wrd_treatment flowsheet,
        Cin=1.2479,
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_inputs_8_19_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)


@pytest.mark.component
def test_TSRO_8_19_21():
    # Stage 3
    # expected_power = 29.3 * pyunits.kW
    expected_power = 20 * pyunits.kW
    m = main(
        Qin=384,
        Cin=4.847 / 2,
        Tin=302,
        Pin=112.6 * pyunits.psi,  # Suction pressure (includes headloss)
        stage_num=3,
        file="wrd_inputs_8_19_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)


@pytest.mark.component
def test_ro_stages_3_13_21():
    # Stage 1
    expected_power = 189.6 * pyunits.kW
    m = main(
        Qin=2451.2,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        stage_num=1,
        file="wrd_inputs_3_13_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)

    # Stage 2
    expected_power = 22.8 * pyunits.kW
    m = m = main(
        Qin=1047.4,
        Cin=1.2479,
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_inputs_3_13_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)


@pytest.mark.component
def test_TSRO_3_13_21():
    # Stage 3
    expected_power = 24.9 * pyunits.kW
    m = main(
        Qin=506.5,
        Cin=4.847 / 2,
        Tin=302,
        Pin=106.3 * pyunits.psi,  # Suction pressure (includes headloss)
        stage_num=3,
        file="wrd_inputs_3_13_21.yaml",
    )
    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
