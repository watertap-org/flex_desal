import pytest
from pyomo.environ import units as pyunits
from wrd.components.ro_system import main
from wrd.utilities import load_config, get_config_value, get_config_file


@pytest.mark.component
def test_ro_system_2_20():
    number_trains = 1
    number_stages = 3
    expected_power = [201.7, 9.7, 22.65]
    expected_perm_flow_gpm = [1467, 648.1, 250.7]
    powers_kW, perm_flows_gpm = main(number_trains, number_stages, date="2_20")
    for t in range(1, number_trains + 1):
        for s in range(1, number_stages + 1):
            # CURRENTLY AVOIDING THIRD STAGE BECAUSE IT DOESN'T MATCH
            modeled_power = powers_kW[f"train_{t}_stage_{s}"]
            assert modeled_power == pytest.approx(
                expected_power[s - 1], rel=0.15
            ), f"Train{t}, Stage {s}: Expected {expected_power[s - 1]} kW, but got {modeled_power} kW"
            modeled_flow = perm_flows_gpm[f"train_{t}_stage_{s}"]
            assert modeled_flow == pytest.approx(
                expected_perm_flow_gpm[s - 1], rel=0.15
            ), f"Train{t}, Stage {s}: Expected {expected_perm_flow_gpm[s - 1]} gpm, but got {modeled_flow} gpm"


# @pytest.mark.component
# def test_ro_system_8_19():
#     number_trains = 1
#     number_stages = 3
#     expected_power = [196.25, 22.71, 29.3]
#     expected_perm_flow_gpm = [1608, 635, 198]
#     powers_kW, perm_flows_gpm = main(number_trains,number_stages,date="8_19")
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
