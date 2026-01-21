import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main


# BECAUSE TOTAL FLOWRATE NOW IN YAML, ONLY THE FULL 4 TRAINS CASE WILL WORK
@pytest.mark.parametrize("num_pro_trains", [4])
@pytest.mark.component
def test_wrd_treatment_train_3_13_21(num_pro_trains):
    file = "wrd_inputs_3_13_21.yaml"
    m = main(
        num_pro_trains=num_pro_trains,
        num_tsro_trains=3,
        uf_split_fraction=[0.4, 0.4, 0.2],
        file=file,
    )

    if num_pro_trains == 4:
        # Expected chemical addition flows are from '2108 GRIP CHEMICALS REPORT.pdf'
        # March 13, 2021

        agg_cost_dict = {
            "sodium_hypochlorite": 9488 * pyunits.USD_2021 / pyunits.month,
            "calcium_hydroxide": 11153 * pyunits.USD_2021 / pyunits.month,
            "sodium_bisulfite": 13835 * pyunits.USD_2021 / pyunits.month,
            "scale_inhibitor": 13568 * pyunits.USD_2021 / pyunits.month,
            "sulfuric_acid": 14376 * pyunits.USD_2021 / pyunits.month,
            "sodium_hydroxide": 29233 * pyunits.USD_2021 / pyunits.month,
            "ammonium_sulfate": 4997 * pyunits.USD_2021 / pyunits.month,
            "electricity": 132498 * pyunits.USD_2021 / pyunits.month,
        }

        for k, v in m.fs.costing.aggregate_flow_costs.items():
            agg_cost = pyunits.convert(v, to_units=pyunits.USD_2021 / pyunits.month)
            # Expected deviations due to daily variations in chemical use
            # i.e., comparing a modeled daily cost to an actual monthly cost
            assert pytest.approx(value(agg_cost), rel=0.35) == value(agg_cost_dict[k])

        expected_LCOW = 0.395 * pyunits.USD_2021 / pyunits.m**3
        assert pytest.approx(value(m.fs.costing.LCOW), rel=0.35) == value(expected_LCOW)

        # Brine
        expected_brine_flow_cost = 61362 * pyunits.USD_2021 / pyunits.month
        assert pytest.approx(
            value(
                pyunits.convert(
                    m.fs.disposal.unit.costing.variable_operating_cost,
                    to_units=pyunits.USD_2021 / pyunits.month,
                )
            ),
            rel=0.35,
        ) == value(expected_brine_flow_cost)
        # Feedwater
        expected_feed_flow_cost = 205820 * pyunits.USD_2021 / pyunits.month
        assert pytest.approx(
            value(
                pyunits.convert(
                    m.fs.feed.costing.variable_operating_cost,
                    to_units=pyunits.USD_2021 / pyunits.month,
                )
            ),
            rel=0.35,
        ) == value(expected_feed_flow_cost)