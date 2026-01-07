import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main


@pytest.mark.unit
def test_no_wrd_file():
    with pytest.raises(
        ValueError, match="Input file must be provided to build WRD system."
    ):
        _ = main(num_pro_trains=2, file=None)


@pytest.mark.parametrize("num_pro_trains", [1, 2, 3, 4])
@pytest.mark.component
def test_wrd_treatment_train_8_19_21(num_pro_trains):
    file = "wrd_inputs_8_19_21.yaml"
    m = main(num_pro_trains=num_pro_trains, file=file)

    if num_pro_trains == 4:
        # Expected chemical addition flows are from '2108 GRIP CHEMICALS REPORT.pdf'
        # for August 19, 2021
        chem_flow_dict = {
            "sodium_hypochlorite": 871 * pyunits.gallon / pyunits.day,
            "calcium_hydroxide": 555 * pyunits.gallon / pyunits.day,
            "sodium_bisulfite": 274 * pyunits.gallon / pyunits.day,
            "scale_inhibitor": 69 * pyunits.gallon / pyunits.day,
            "sulfuric_acid": 365 * pyunits.gallon / pyunits.day,
            "sodium_hydroxide": 304 * pyunits.gallon / pyunits.day,
            "ammonium_sulfate": 130 * pyunits.gallon / pyunits.day,
        }
        for k, v in chem_flow_dict.items():
            blk = m.fs.find_component(f"{k}_addition")
            chem_flow = pyunits.convert(
                blk.unit.chemical_soln_flow_vol, to_units=pyunits.gallon / pyunits.day
            )
            assert pytest.approx(value(chem_flow), rel=0.15) == value(v)

        # Expected chemical costs are from ARC Cost Tracker, 'Chem' tab, entry for August 2021
        # Expected ekectricity cost from Arc Cost Tracker, 'SCE' tab, entry for August 2021 (column R)
        agg_cost_dict = {
            "sodium_hypochlorite": 17088 * pyunits.USD_2021 / pyunits.month,
            "calcium_hydroxide": 16939 * pyunits.USD_2021 / pyunits.month,
            "sodium_bisulfite": 20940 * pyunits.USD_2021 / pyunits.month,
            "scale_inhibitor": 14606 * pyunits.USD_2021 / pyunits.month,
            "sulfuric_acid": 15086 * pyunits.USD_2021 / pyunits.month,
            "sodium_hydroxide": 7814 * pyunits.USD_2021 / pyunits.month,
            "ammonium_sulfate": 6061 * pyunits.USD_2021 / pyunits.month,
            "electricity": 197223 * pyunits.USD_2021 / pyunits.month,
        }
        for k, v in m.fs.costing.aggregate_flow_costs.items():
            agg_cost = pyunits.convert(v, to_units=pyunits.USD_2021 / pyunits.month)
            # Expected deviations due to daily variations in chemical use
            # i.e., comparing a modeled daily cost to an actual monthly cost
            assert pytest.approx(value(agg_cost), rel=0.35) == value(agg_cost_dict[k])

        expected_LCOW = 0.442 * pyunits.USD_2021 / pyunits.m**3
        assert pytest.approx(value(m.fs.costing.LCOW), rel=0.35) == value(expected_LCOW)

    # TODO: Add more tests against facility data


@pytest.mark.parametrize("num_pro_trains", [1, 2, 3, 4])
@pytest.mark.component
def test_wrd_treatment_train_3_13_21(num_pro_trains):
    file = "wrd_inputs_3_13_21.yaml"
    m = main(num_pro_trains=num_pro_trains, file=file)
    # TODO: Add tests against facility data
