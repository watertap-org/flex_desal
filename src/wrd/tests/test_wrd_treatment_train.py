import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main
from wrd.utilities import load_config, get_config_value, get_config_file


@pytest.mark.parametrize("num_pro_trains", [4, 3, 2])
@pytest.mark.component
def test_wrd_treatment_train(num_pro_trains):
    m = main(num_pro_trains=num_pro_trains)

@pytest.fixture(scope="module")
def full_wrd_system_8_19_21():
    m = main(num_pro_trains=1, num_tsro_trains=1, num_pro_stages=2)
    return m


def _get_stage_objects(m, train_idx, stage_idx):
    if stage_idx < 3:
        train = m.fs.train[train_idx]
        stage = train.stage[stage_idx]
    else:
        stage = m.fs.tsro_train[train_idx]
    pump = stage.pump
    perm_flow = stage.ro.unit.mixed_permeate[0].flow_vol_phase["Liq"]
    return pump, perm_flow


# @pytest.mark.component
# def test_with_1_train():
#     m = main(num_pro_trains=1, num_tsro_trains=1, num_pro_stages=2)


# @pytest.mark.component
# def test_with_4_train():
#     m = main(num_pro_trains=4, num_tsro_trains=4, num_pro_stages=2)
    return m


def _get_stage_objects(m, train_idx, stage_idx):
    if stage_idx < 3:
        train = m.fs.train[train_idx]
        stage = train.stage[stage_idx]
    else:
        stage = m.fs.tsro_train[train_idx]
    pump = stage.pump
    perm_flow = stage.ro.unit.mixed_permeate[0].flow_vol_phase["Liq"]
    return pump, perm_flow


# @pytest.mark.component
# def test_with_1_train():
#     m = main(num_pro_trains=1, num_tsro_trains=1, num_pro_stages=2)


# @pytest.mark.component
# def test_with_4_train():
#     m = main(num_pro_trains=4, num_tsro_trains=4, num_pro_stages=2)


# Current test only checks first stage of first train pump power and perm flow
# Current test only checks first stage of first train pump power and perm flow
@pytest.mark.component
def test_wrd_treatment_train_PRO_pump_1(full_wrd_system_8_19_21):
    pump, _ = _get_stage_objects(full_wrd_system_8_19_21, 1, 1)
    power = pyunits.convert(pump.unit.work_mechanical[0], to_units=pyunits.kW)
    power = pyunits.convert(pump.unit.work_mechanical[0], to_units=pyunits.kW)
    expected_power = 196.25 * pyunits.kW
    # Units check
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)


# Current test only checks first strage of first train pump power
@pytest.mark.component
def test_wrd_treatment_train_PRO1(full_wrd_system_8_19_21):
    _, perm_flow = _get_stage_objects(full_wrd_system_8_19_21, 1, 1)
    expected_perm_flow = 1608.2 * pyunits.gal / pyunits.minute
    expected_perm_flow = pyunits.convert(
        expected_perm_flow, to_units=pyunits.m**3 / pyunits.s
    )
    # Units check
    assert_units_consistent(perm_flow + expected_perm_flow)
    assert value(perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


# Will eventually add more tests for other pumps and total system power, etc.
@pytest.mark.skip
def test_wrd_treatment_train_total_power(full_wrd_system_8_19_21):
    pump, _ = _get_stage_objects(full_wrd_system_8_19_21, 1, 1)
    power = pyunits.convert(pump.pump.work_mechanical[0], to_units=pyunits.kW)
    expected_power = 196.25 * pyunits.kW
    # Units check
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)
