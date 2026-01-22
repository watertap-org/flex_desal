import pytest

from pyomo.environ import value, units as pyunits

import wrd.components.ro_stage as ro_stage
import wrd.wrd_treatment_train as wrd_full_sys
import wrd.components.ro_train as ro_train

# RO stage tests
@pytest.mark.skip
def test_ro_PRO1_8_19_21():
    # Stage 1
    m = ro_stage.main(
        Qin=2637, # <--- measured value
        Cin=0.528, # <--- measured value
        Tin=302, # Why are we testing everything at this temp and not 298 K?
        Pin=35.4 * pyunits.psi, # <--- measured value
        stage_num=1,
        file="wrd_inputs_8_19_21.yaml",
    )

    expected_power = 196.25 * pyunits.kW  # <--- measured value
    expected_perm_flow = 1608.2 * pyunits.gal / pyunits.min # <--- measured value

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)


@pytest.mark.skip
def test_ro_PRO2_8_19_21():
    m = ro_stage.main(
        Qin=1029,  
        Cin=1.2479, 
        Tin=302,
        Pin=131.2 * pyunits.psi,
        stage_num=2,
        file="wrd_inputs_8_19_21.yaml",
    )

    expected_power = 22.71 * pyunits.kW  # <--- measured value
    expected_perm_flow = 635 * pyunits.gal / pyunits.min # <--- measured value

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )

    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)
    # Add permeate salinity check?

@pytest.mark.skip
def test_TSRO_8_19_21():
    m = ro_stage.main(
        Qin=384,
        Cin= 2.4235,
        Tin=302,
        Pin=112.6 * pyunits.psi,  # Suction pressure (includes headloss)
        stage_num=3,
        file="wrd_inputs_8_19_21.yaml",
    )

    expected_power = 29.3 * pyunits.kW  # <--- measured value
    expected_perm_flow = 198 * pyunits.gal / pyunits.min # <--- measured value

    actual_power = pyunits.convert(
        m.fs.ro_stage.pump.unit.work_mechanical[0], to_units=pyunits.kW
    )
    actual_perm_flow = pyunits.convert(
        m.fs.ro_stage.ro.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )
    assert value(actual_power) == pytest.approx(value(expected_power), rel=0.15)
    assert value(actual_perm_flow) == pytest.approx(value(expected_perm_flow), rel=0.15)

# Arguably these tests are overkill because they are just testing the stages again
@pytest.mark.component
def test_ro_train1_8_19_21():
    expected_power = (196.25 + 22.7) * pyunits.kW
    expected_product_flow = (1608 + 635) * pyunits.gal / pyunits.min
    expected_SEC = pyunits.convert(
        expected_power / expected_product_flow, to_units=pyunits.kWh / pyunits.m**3
    )
    m = ro_train.main()

    actual_power = pyunits.convert(m.fs.ro_train.total_pump_power, to_units=pyunits.kW)
    assert pytest.approx(value(actual_power), rel=0.15) == value(expected_power)

    actual_product_flow = pyunits.convert(
        m.fs.ro_train.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.gal / pyunits.min,
    )
    assert pytest.approx(value(actual_product_flow), rel=0.15) == value(
        expected_product_flow
    )
    assert pytest.approx(value(m.fs.costing.SEC), rel=0.15) == value(expected_SEC)


# BECAUSE TOTAL FLOWRATE NOW IN YAML, ONLY THE FULL 4 TRAINS CASE WILL WORK
@pytest.mark.parametrize("num_pro_trains", [4])
@pytest.mark.skip
def test_wrd_treatment_train_8_19_21(num_pro_trains):
    file = "wrd_inputs_8_19_21.yaml"
    m = wrd_full_sys.main(num_pro_trains=num_pro_trains, uf_split_fraction=[0.4, 0.4, 0.2], file=file)

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

    # Brine
    expected_brine_flow_cost = 59162 * pyunits.USD_2021 / pyunits.month
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
    expected_feed_flow_cost = 241400 * pyunits.USD_2021 / pyunits.month
    assert pytest.approx(
        value(
            pyunits.convert(
                m.fs.feed.costing.variable_operating_cost,
                to_units=pyunits.USD_2021 / pyunits.month,
            )
        ),
        rel=0.35,
    ) == value(expected_feed_flow_cost)
    # TODO: Add more tests against facility data.
    # power of some pumps, or total plant power use.

