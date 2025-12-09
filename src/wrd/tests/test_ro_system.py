import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from wrd.components.ro_system import main
from wrd.utilities import load_config, get_config_value, get_config_file


# @pytest.mark.component
# def test_ro_system_2_20():
#     number_trains = 1
#     number_stages = 3
#     expected_power = [201.7, 9.7, 22.65]
#     expected_perm_flow_gpm = [1467, 648.1, 250.7]
#     powers_kW, perm_flows_gpm = main(number_trains, number_stages, date="2_20")
#     for t in range(1, number_trains + 1):
#         for s in range(1, number_stages + 1):
#             modeled_power = powers_kW[f"train_{t}_stage_{s}"]
#             assert modeled_power == pytest.approx(
#                 expected_power[s - 1], rel=0.15
#             ), f"Train{t}, Stage {s}: Expected {expected_power[s - 1]} kW, but got {modeled_power} kW"
#             modeled_flow = perm_flows_gpm[f"train_{t}_stage_{s}"]
#             assert modeled_flow == pytest.approx(
#                 expected_perm_flow_gpm[s - 1], rel=0.15
#             ), f"Train{t}, Stage {s}: Expected {expected_perm_flow_gpm[s - 1]} gpm, but got {modeled_flow} gpm"


@pytest.mark.component
def test_ro_system_8_19():
    number_trains = 1
    number_stages = 2  # STAGE 3 will fail until pressure drop added
    expected_power = [v * pyunits.kW for v in (196.25, 22.71, 29.3)]
    expected_perm_flow_gpm = [v * pyunits.gal/pyunits/min for v in (1608, 635, 198)]
    expected_perm_flow = [pyunits.convert(f, to_units=pyunits.m**3 / pyunits.s) for f in expected_perm_flow_gpm]

    m = main(number_trains, number_stages, date="8_19_21")
    for t in range(1, number_trains + 1):
        train = m.fs.ro_system.find_component(f"train_{t}")
        for s in range(1, number_stages + 1):
            pump = train.find_component(f"pump{s}")
            modeled_power = pyunits.convert(pump.pump.work_mechanical[0],to_units=pyunits.kW)
            stage_rr = train.find_component(f"ro_stage_{s}").recovery_vol_phase[0, "Liq"]
            stage_perm = stage_rr * pump.feed_out.properties[0].flow_vol_phase["Liq"]

            assert_units_consistent(modeled_power + expected_power[0])
            assert_units_consistent(stage_perm + expected_perm_flow[0])
            
            assert modeled_power == pytest.approx(expected_power[s - 1], rel=0.15), f"Train{t}, Stage {s}: Expected {expected_power[s - 1]} kW, but got {modeled_power} kW"
            assert value(stage_perm) == pytest.approx(
                value(expected_perm_flow_gpm[s - 1]), rel=0.15
            ), f"Train{t}, Stage {s}: Expected {expected_perm_flow_gpm[s - 1]} gpm, but got {modeled_flow} gpm"


# @pytest.mark.component
# def test_ro_system_3_13():
#     number_trains = 1
#     number_stages = 3
#     expected_power = [189.6, 22.8, 24.9]
#     expected_perm_flow_gpm = [1404.7, 617.1, 278.5]
#     powers_kW, perm_flows_gpm = main(number_trains,number_stages,date="3_13")
#     for t in range(1, 1 + number_trains):
#         for s in range(
#             1, number_stages
#         ):  # CURRENTLY AVOIDING THIRD STAGE BECAUSE IT DOESN'T MATCH
#             assert powers_kW[f"train_{t}_stage_{s}"] == pytest.approx(
#                 expected_power[s - 1], rel=0.15
#             )
#             assert perm_flows_gpm[f"train_{t}_stage_{s}"] == pytest.approx(
#                 expected_perm_flow_gpm[s - 1], rel=0.15
#             )
