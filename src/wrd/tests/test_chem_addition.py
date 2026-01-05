import pytest

from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_equivalent

import wrd.components.chemical_addition as ca


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
    # These flows are based on Qin=10000 and yaml inputs (see Chemical_cost_data.xlsx in Box)
    vol_flow_rates = {
        "ammonium_sulfate": 6.41e-7,
        "sodium_hypochlorite": 1.69e-5,
        "sulfuric_acid": 3.55e-5,
        "scale_inhibitor": 2.58e-6,
        "calcium_hydroxide": 6.52e-5,
        "sodium_hydroxide": 5.96e-6,
        "sodium_bisulfite": 4.78e-6,
    }
    for chem, flow in vol_flow_rates.items():
        m = ca.main(
            chemical_name=chem,
            Qin=10000,
            Cin=0.5,
            dose=None,
            chem_cost=None,
        )
        chem_vol_flow = m.fs.chem_addition.unit.chemical_soln_flow_vol
        expected_vol_flow = flow * pyunits.m**3 / pyunits.s
        assert_units_equivalent(chem_vol_flow, expected_vol_flow)
        assert pytest.approx(value(chem_vol_flow), rel=0.05) == value(
            expected_vol_flow
        )  # m^3/s


@pytest.mark.component
def test_chem_addition_costs():
    # These costs are based directly on Qin=10000 and yaml inputs, so values should agree very closely
    monthly_costs = {
        "ammonium_sulfate": 841,
        "sodium_hypochlorite": 8603,
        "sulfuric_acid": 35611,
        "scale_inhibitor": 15626,
        "calcium_hydroxide": 143327,
        "sodium_hydroxide": 11810,
        "sodium_bisulfite": 9322,
    }
    for chem, cost in monthly_costs.items():
        m = ca.main(
            chemical_name=chem,
            Qin=10000,
            Cin=0.5,
            dose=None,
            chem_cost=None,
        )
        operational_cost = m.fs.costing.aggregate_flow_costs[chem]
        assert_units_equivalent(
            operational_cost, m.fs.costing.base_currency / m.fs.costing.base_period
        )
        assert pytest.approx(value(operational_cost), rel=0.05) == cost  # $/month
