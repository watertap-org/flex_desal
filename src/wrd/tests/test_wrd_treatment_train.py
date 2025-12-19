import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main
from pyomo.util.check_units import assert_units_consistent


# Parametrized fixture for model creation
@pytest.fixture(params=[4, 3, 2], scope="module")
def wrd_treatment_train_model(request):
    num_pro_trains = request.param
    m = main(num_pro_trains=num_pro_trains)
    return m, num_pro_trains


@pytest.mark.component
def test_wrd_treatment_train(wrd_treatment_train_model):
    m, num_pro_trains = wrd_treatment_train_model
    # Optionally, add assertions here to check model properties


def _get_stage_objects(m, train_idx, stage_idx):
    if stage_idx < 3:
        train = m.fs.train[train_idx]
        stage = train.stage[stage_idx]
    else:
        stage = m.fs.tsro_train[train_idx]
    pump = stage.pump
    perm_flow = stage.ro.unit.mixed_permeate[0].flow_vol_phase["Liq"]
    return pump, perm_flow


# Current test only checks first stage of first train pump power and perm flow
@pytest.mark.component
def test_wrd_treatment_train_PRO_pump_1(wrd_treatment_train_model):
    m, _ = wrd_treatment_train_model
    pump, _ = _get_stage_objects(m, 1, 1)
    power = pyunits.convert(pump.unit.work_mechanical[0], to_units=pyunits.kW)
    expected_power = 196.25 * pyunits.kW
    # Units check
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)


# Current test only checks first stage of first train pump power
@pytest.mark.component
def test_wrd_treatment_train_PRO1(wrd_treatment_train_model):
    m, _ = wrd_treatment_train_model
    _, perm_flow = _get_stage_objects(m, 1, 1)
    expected_perm_flow = 1608.2 * pyunits.gal / pyunits.minute
    expected_perm_flow = pyunits.convert(
        expected_perm_flow, to_units=pyunits.m**3 / pyunits.s
    )
    # Units check
    assert_units_consistent(perm_flow + expected_perm_flow)
    assert value(perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


# Not currently working as total system power hasn't been cacluated from data.
@pytest.mark.skip
def test_wrd_treatment_train_total_power(wrd_treatment_train_model):
    m, _ = wrd_treatment_train_model
    power = pyunits.convert(m.fs.total_system_pump_power, to_units=pyunits.kW)
    expected_power = 196.25 * pyunits.kW
    # Units check
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)
