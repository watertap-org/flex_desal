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
    mass_flow_rates = {
        "ammonium_sulfate": 8.32e-5,
        "sodium_hypochlorite": 6.65e-4,
        "sulfuric_acid": 1.60e-2,
        "scale_inhibitor": 7.49e-4,
        "calcium_hydroxide": 7.82e-3,
        "sodium_hydroxide": 4.99e-4,
        "sodium_bisulfite": 6.65e-4,
    }
    for chem, flow in mass_flow_rates.items():
        m = ca.main(
            chemical_name=chem,
            Qin=2637,
            Cin=0.5,
            dose=None,
            chem_cost=None,
            chem_purity=None,
        )

        chem_mass_flow = m.fs.costing.find_component(f"aggregate_flow_{chem}")
        expected_mass_flow = flow * pyunits.kg / pyunits.s
        assert pytest.approx(value(chem_mass_flow), rel=0.05) == value(
            expected_mass_flow
        )  # kg/s


@pytest.mark.component
def test_chem_addition_costs():
    # These costs are based directly on Qin=2637 and yaml inputs, so values should agree very closely
    monthly_costs = {
        "ammonium_sulfate": 222,
        "sodium_hypochlorite": 2256,
        "sulfuric_acid": 8880,
        "scale_inhibitor": 4533,
        "calcium_hydroxide": 47309,
        "sodium_hydroxide": 3114,
        "sodium_bisulfite": 2458,
    }
    for chem, cost in monthly_costs.items():
        m = ca.main(
            chemical_name=chem,
            Qin=2637,
            Cin=0.5,
            dose=None,
            chem_cost=None,
            chem_purity=None,
        )

        operational_cost = m.fs.costing.aggregate_flow_costs[chem]
        assert_units_equivalent(
            operational_cost, m.fs.costing.base_currency / m.fs.costing.base_period
        )
        assert pytest.approx(value(operational_cost), rel=0.05) == cost  # $/month
