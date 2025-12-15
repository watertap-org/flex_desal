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

def _get_chem_cost(m,chem_name):
    chem_cost = m.fs.costing.find_component(f"chemical_cost_{chem_name}").cost
    return chem_cost

# Current test only checks first strage of first train pump power
@pytest.mark.component
def test_wrd_treatment_train_PRO_pump_1(full_wrd_system_8_19_21):
    pump, _ = _get_stage_objects(full_wrd_system_8_19_21, 1, 1)
    power = pyunits.convert(pump.pump.work_mechanical[0], to_units=pyunits.kW)
    expected_power = 196.25 * pyunits.kW
    # Units check
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)

@pytest.mark.component
def test_wrd_treatment_ammonium_sulfate(full_wrd_system_8_19_21):
    chem_cost = _get_chem_cost(full_wrd_system_8_19_21, "ammonium_sulfate")
    expected_chemical_cost = 1 * pyunits.USD_2018
    # Units check
    assert_units_consistent(chem_cost + expected_chemical_cost)
    assert value(chem_cost) == pytest.approx(value(expected_chemical_cost), rel=0.15)

# Can copy and paste the above for other specific chemicals / pumps / components


# Total Costs

@pytest.mark.component
def test_wrd_treatment_train_total_power(full_wrd_system_8_19_21):
    power = pyunits.convert(full_wrd_system_8_19_21.fs.total_power, to_units=pyunits.kW)
    expected_power = 1 * pyunits.kW
    # Units check
    assert_units_consistent(power + expected_power)
    assert value(power) == pytest.approx(value(expected_power), rel=0.15)

@pytest.mark.component
def test_wrd_treatment_chemical_costs(full_wrd_system_8_19_21):
    expected_chemical_cost = 1 * pyunits.USD_2018
    # Units check
    assert_units_consistent(full_wrd_system_8_19_21.fs.chemical_costs + expected_chemical_cost)
    assert value(full_wrd_system_8_19_21.fs.chemical_costs) == pytest.approx(value(expected_chemical_cost), rel=0.15)


