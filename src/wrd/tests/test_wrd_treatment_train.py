import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from wrd.wrd_treatment_train import main
from wrd.utilities import load_config, get_config_value, get_config_file


@pytest.fixture(scope="module")
def full_wrd_system_8_19_21():
    number_stages = 3
    m = main(number_stages=number_stages, date="8_19_21")
    return m


def _get_stage_objects(m, train_idx, stage_idx):
    train = m.fs.ro_system.find_component(f"train_{train_idx}")
    pump = train.find_component(f"pump{stage_idx}")
    stage = train.find_component(f"ro_stage_{stage_idx}")
    perm_flow = stage.mixed_permeate[0].flow_vol_phase["Liq"]
    return pump, perm_flow


# @pytest.mark.component
# def test_wrd_treatment_train_builds_8_19_21(full_wrd_system_8_19_21):
#     pump, _ = _get_stage_objects(full_wrd_system_8_19_21, 1, 1)
#     power = pyunits.convert(pump.pump.work_mechanical[0], to_units=pyunits.kW)
#     expected_power = 196.25 * pyunits.kW
#     # Units check
#     assert_units_consistent(power + expected_power)
#     assert value(power) == pytest.approx(value(expected_power), rel=0.15)
