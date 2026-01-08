import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.ro_stage import main


@pytest.mark.component
def test_ro_PRO1_8_19_21():
    # Stage 1
    m = main(
        Qin=2637,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        stage_num=1,
        file="wrd_inputs_8_19_21.yaml",
    )

    expected_power = 196.25 * pyunits.kW
    expected_perm_flow = 1608.2 * pyunits.gal / pyunits.min

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


@pytest.mark.component
def test_ro_PRO2_8_19_21():
    m = main(
        Qin=1029,  # -> measured value, 1175 gpm is flow in modeled in wrd_treatment flowsheet
        Cin=1.2479,
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_inputs_8_19_21.yaml",
    )

    # expected_power = 22.71 * pyunits.kW <--- measured value
    expected_power = 20 * pyunits.kW  # <--- modeled value
    expected_perm_flow = 635 * pyunits.gal / pyunits.min

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


@pytest.mark.component
def test_TSRO_8_19_21():
    m = main(
        Qin=384,
        Cin=4.847 / 2,
        Tin=302,
        Pin=112.6 * pyunits.psi,  # Suction pressure (includes headloss)
        stage_num=3,
        file="wrd_inputs_8_19_21.yaml",
    )

    # expected_power = 29.3 * pyunits.kW <--- measured value
    expected_power = 20 * pyunits.kW  # <--- modeled value
    # expected_perm_flow = 198 * pyunits.gal / pyunits.min <--- measured value
    expected_perm_flow = 267 * pyunits.gal / pyunits.min  # <--- modeled value

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


@pytest.mark.component
def test_ro_PRO1_3_13_21():
    # Stage 1
    m = main(
        Qin=2451.2,
        Cin=0.528,
        Tin=302,
        Pin=101325,
        stage_num=1,
        file="wrd_inputs_3_13_21.yaml",
    )

    expected_power = 189.6 * pyunits.kW
    # expected_perm_flow = 1404.7 * pyunits.gal / pyunits.min # measured value
    expected_perm_flow = 1797 * pyunits.gal / pyunits.min  # modeled value

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


@pytest.mark.component
def test_ro_PRO2_3_13_21():
    # Stage 2
    m = main(
        Qin=1047.4,
        Cin=1.2479,
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_inputs_3_13_21.yaml",
    )

    expected_power = 22.8 * pyunits.kW
    expected_perm_flow = 617.1 * pyunits.gal / pyunits.min

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


@pytest.mark.component
def test_TSRO_3_13_21():
    # Stage 3
    m = main(
        Qin=506.5,
        Cin=4.847 / 2,
        Tin=302,
        Pin=106.3 * pyunits.psi,  # Suction pressure (includes headloss)
        stage_num=3,
        file="wrd_inputs_3_13_21.yaml",
    )

    # expected_power = 24.9 * pyunits.kW <--- measured value
    expected_power = 20 * pyunits.kW  # <--- modeled value
    expected_perm_flow = 278.5 * pyunits.gal / pyunits.min

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)
