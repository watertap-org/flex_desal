import pytest

from pyomo.environ import assert_optimal_termination, value, units as pyunits

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError

from watertap.core.solvers import get_solver

import wrd.components.chemical_addition as ca

solver = get_solver()


@pytest.mark.component
def test_chem_addition_chem_flow():
    # These costs are based directly on Qin=2637 and yaml inputs, so values should agree very closely
    mass_flow_rates = {
        "ammonium_sulfate": 8.32e-5,
        "sodium_hypochlorite": 6.65e-4,
        "sulfuric_acid": 1.60e-2,
        "scale_inhibitor": 7.49e-4,
        "calcium_hydroxide": 7.82e-3,
        "sodium_hydroxide": 4.99e-4,
        "sodium_bisulfite": 6.65e-4,
    }
    for i, chem in enumerate(mass_flow_rates.keys(), 1):
        m = ca.main(
            chemical_name=chem,
            Qin=2637,
            Cin=0.5,
            dose=None,
            chem_cost=None,
            chem_purity=None,
        )

        chem_mass_flow = m.fs.costing.find_component(f"aggregate_flow_{chem}")
        expected_mass_flow = mass_flow_rates[chem] * pyunits.kg / pyunits.s
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
    for i, chem in enumerate(monthly_costs.keys(), 1):
        # All of this stuff should be in main??

        # Dummy data just for testing
        # dose = None, #0.01 * i
        # cost = None, #0.5 * i
        # purity = None, #1

        # m = ca.build_system(chemical_name=chem)
        # ca.calculate_scaling_factors(m)
        # ca.set_inlet_conditions(m)
        # ca.set_chem_addition_op_conditions(m.fs.chem_addition, dose=dose)
        # ca.initialize_system(m)
        # assert degrees_of_freedom(m) == 0
        # results = solver.solve(m)
        # assert_optimal_termination(results)
        # ca.report_chem_addition(m.fs.chem_addition, w=40)
        # chem_registered = chem in m.fs.costing._registered_flows.keys()
        # ca.add_chem_addition_costing(
        #     m.fs.chem_addition, chem_cost=cost, chem_purity=purity
        # )
        # assert not chem_registered
        # m.fs.costing.cost_process()
        # m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        # m.fs.costing.initialize()
        # assert degrees_of_freedom(m) == 0
        # results = solver.solve(m)
        # assert_optimal_termination(results)

        m = ca.main(
            chemical_name=chem,
            Qin=2637,
            Cin=0.5,
            dose=None,
            chem_cost=None,
            chem_purity=None,
        )

        operational_cost = value(m.fs.costing.aggregate_flow_costs[chem])
        expected_cost = monthly_costs[chem]
        assert pytest.approx(value(operational_cost), rel=0.05) == expected_cost  # $/month


# Don't understand what this is testing-->
@pytest.mark.skip
def test_chem_addition_missing_data():
    msg = "Must specify a chemical for addition."
    with pytest.raises(ConfigurationError, match=msg):
        m = ca.build_system(chemical_name=None)
