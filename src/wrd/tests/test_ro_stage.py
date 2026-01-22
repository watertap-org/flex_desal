import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.ro_stage import main


@pytest.mark.component
def test_ro_PRO_stage_1():
    # Stage 1
    m = main(
        Qin=2451.2,
        Cin=0.528,
        Tin=302,
        Pin=35.4 * pyunits.psi,
        stage_num=1,
        file="wrd_inputs_3_13_21.yaml",
    )

    expected_power = 180.74 * pyunits.kW
    expected_perm_flow = 1413.79 * pyunits.gal / pyunits.min  # modeled value

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=1e-3)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=1e-3)


@pytest.mark.component
def test_ro_PRO_stage_2():
    # Stage 2
    m = main(
        Qin=1047.4,
        Cin=1.2479,
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_inputs_3_13_21.yaml",
    )

    expected_power = 25.97 * pyunits.kW
    expected_perm_flow = 625.26 * pyunits.gal / pyunits.min

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=1e-3)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=1e-3)


@pytest.mark.component
def test_TSRO():
    # Stage 3
    m = main(
        Qin=506.5,
        Cin=4.847 / 2,
        Tin=302,
        Pin=106.3 * pyunits.psi,  # Suction pressure (includes headloss)
        stage_num=3,
        file="wrd_inputs_3_13_21.yaml",
    )

    expected_power = 19.43 * pyunits.kW
    expected_perm_flow = 287.21 * pyunits.gal / pyunits.min

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=1e-3)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=1e-3)
    # Add permeate salinity check?