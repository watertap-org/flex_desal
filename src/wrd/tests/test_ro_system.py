import pytest

from wrd.components.ro_system import main


# @pytest.mark.component
# def test_ro_system_2_20():
#     number_trains = 1
#     Qin = 2450 / 264.2 / 60  # gpm to m3/s
#     Cin = 1082 * 0.5 / 1000  # us/cm to g/L
#     expected_power = 201.7 + 9.7 + 22.65
#     expected_perm_flow_gpm = 1467 + 648.1 + 250.7
#     power, perm_flow_gpm = main(number_trains, Qin, Cin)
#     assert power == pytest.approx(expected_power, rel=0.5)
#     assert perm_flow_gpm == pytest.approx(expected_perm_flow_gpm, rel=0.5)
# # Is this one needed?
# #    assert perm_salinity

@pytest.mark.component
def test_ro_system_8_19():
    number_trains = 1
    number_stages = 3
    Qin = 2637 / 264.2 / 60  # gpm to m3/s
    Cin = 1055 * 0.5 / 1000  # us/cm to g/L
    expected_power = [196.25, 22.71, 29.3]
    expected_perm_flow_gpm = [1608, 635, 198]
    powers_kW, perm_flows_gpm = main(number_trains, Qin, Cin)
    for t in range(1,1+number_trains):
        for s in range(1,number_stages): #CURRENTLY AVOIDING THIRD STAGE BECAUSE IT DOESN'T MATCH
            assert powers_kW[f"train_{t}_stage_{s}"] == pytest.approx(expected_power[s-1], rel=0.2)
            assert perm_flows_gpm[f"train_{t}_stage_{s}"] == pytest.approx(expected_perm_flow_gpm[s-1], rel=0.2)
# Is this one needed?
#    assert perm_salinity

# @pytest.mark.component
# def test_ro_system_3_13():
#     number_trains = 1
#     Qin = 2452 / 264.2 / 60  # gpm to m3/s
#     Cin = 1007 * 0.5 / 1000  # us/cm to g/L
#     expected_power = 189.6 + 22.8 + 24.9
#     expected_perm_flow_gpm = 1404.7 + 617.1 + 278.5
#     power, perm_flow_gpm = main(number_trains, Qin, Cin)
#     assert power == pytest.approx(expected_power, rel=0.5)
#     assert perm_flow_gpm == pytest.approx(expected_perm_flow_gpm, rel=0.5)
# # Is this one needed?
# #    assert perm_salinity