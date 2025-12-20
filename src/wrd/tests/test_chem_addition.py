import pytest

from pyomo.environ import assert_optimal_termination

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError

from watertap.core.solvers import get_solver

import wrd.components.chemical_addition as ca

solver = get_solver()


@pytest.mark.component
def test_chem_addition():
    for i, chem in enumerate(
        [
            "ammonium_sulfate",
            "sodium_hypochlorite",
            "sulfuric_acid",
            "scale_inhibitor",
            "calcium_hydroxide",
            "sodium_bisulfite",
        ],
        1,
    ):
        # Dummy data just for testing
        dose = 0.01 * i
        cost = 0.5 * i
        purity = 1

        m = ca.build_system(chemical_name=chem)
        ca.calculate_scaling_factors(m)
        ca.set_inlet_conditions(m)
        ca.set_chem_addition_op_conditions(m.fs.chem_addition, dose=dose)
        ca.initialize_system(m)

        assert degrees_of_freedom(m) == 0
        results = solver.solve(m)
        assert_optimal_termination(results)
        ca.report_chem_addition(m.fs.chem_addition, w=40)

        chem_registered = chem in m.fs.costing._registered_flows.keys()

        ca.add_chem_addition_costing(
            m.fs.chem_addition, chem_cost=cost, chem_purity=purity
        )

        assert not chem_registered

        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0
        results = solver.solve(m)
        assert_optimal_termination(results)


@pytest.mark.component
def test_chem_addition_missing_data():
    msg = "Must specify a chemical for addition."
    with pytest.raises(ConfigurationError, match=msg):
        m = ca.build_system(chemical_name=None)
