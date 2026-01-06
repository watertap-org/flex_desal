import pytest

from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_equivalent

import wrd.components.chemical_addition as ca

flow_pre_treat = 10652
flow_post_treat = 9397


@pytest.mark.unit
def test_chem_addition_missing_data():
    msg = "Must specify a chemical for addition."
    with pytest.raises(ValueError, match=msg):
        _ = ca.build_system(chemical_name=None)


@pytest.mark.component
def test_chem_addition_chem_flow():
    # These flows are based directly on Qin=2637 and yaml inputs, so values should agree very closely
    # mass_flow_rates = {
    #     "ammonium_sulfate": 8.32e-5,
    #     "sodium_hypochlorite": 6.65e-4,
    #     "sulfuric_acid": 1.60e-2,
    #     "scale_inhibitor": 1,
    #     "calcium_hydroxide": 7.82e-3,
    #     "sodium_hydroxide": 4.99e-4,
    #     "sodium_bisulfite": 6.65e-4,
    # }
    # These flows are based on inlet flow for August 19, 2021 WRD treatment train (Qin=10652 gpm)

    chem_flow_dict = {
        "sodium_hypochlorite": 871 * pyunits.gallon / pyunits.day,
        "calcium_hydroxide": 555 * pyunits.gallon / pyunits.day,
        "sodium_bisulfite": 274 * pyunits.gallon / pyunits.day,
        "scale_inhibitor": 69 * pyunits.gallon / pyunits.day,
        "sulfuric_acid": 365 * pyunits.gallon / pyunits.day,
        "sodium_hydroxide": 304 * pyunits.gallon / pyunits.day,
        "ammonium_sulfate": 130 * pyunits.gallon / pyunits.day,
    }
    for chem, v in chem_flow_dict.items():
        if chem in [
            "ammonium_sulfate",
            "sodium_hypochlorite",
            "sulfuric_acid",
            "scale_inhibitor",
        ]:
            flow = flow_pre_treat
        else:
            flow = flow_post_treat
        m = ca.main(
            chemical_name=chem,
            Qin=flow,
        )
        chem_vol_flow = m.fs.chem_addition.unit.chemical_soln_flow_vol
        chem_flow = pyunits.convert(
            chem_vol_flow, to_units=pyunits.gallon / pyunits.day
        )
        assert pytest.approx(value(chem_flow), rel=0.15) == value(v)


@pytest.mark.component
def test_chem_addition_costs():
    # These costs are from August 19, 2021 WRD treatment train
    # This is facility data
    monthly_costs = {
        "sodium_hypochlorite": 17088 * pyunits.USD_2021 / pyunits.month,
        "calcium_hydroxide": 16939 * pyunits.USD_2021 / pyunits.month,
        "sodium_bisulfite": 20940 * pyunits.USD_2021 / pyunits.month,
        "scale_inhibitor": 14606 * pyunits.USD_2021 / pyunits.month,
        "sulfuric_acid": 15086 * pyunits.USD_2021 / pyunits.month,
        "sodium_hydroxide": 6688 * pyunits.USD_2021 / pyunits.month,
        "ammonium_sulfate": 6061 * pyunits.USD_2021 / pyunits.month,
    }
    for chem, cost in monthly_costs.items():
        if chem in [
            "ammonium_sulfate",
            "sodium_hypochlorite",
            "sulfuric_acid",
            "scale_inhibitor",
        ]:
            flow = flow_pre_treat
        else:
            flow = flow_post_treat
        m = ca.main(
            chemical_name=chem,
            Qin=flow,
        )
        operational_cost = m.fs.costing.aggregate_flow_costs[chem]
        agg_cost = pyunits.convert(
            operational_cost, to_units=pyunits.USD_2021 / pyunits.month
        )
        assert_units_equivalent(
            operational_cost, m.fs.costing.base_currency / m.fs.costing.base_period
        )
        assert pytest.approx(value(agg_cost), rel=0.15) == value(cost)  # $/month
