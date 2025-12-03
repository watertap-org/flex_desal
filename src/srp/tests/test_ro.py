import pytest

from pyomo.environ import value, assert_optimal_termination, units as pyunits

from watertap.core.solvers import get_solver

import srp.components.ro as ro


@pytest.mark.component
def test_ro_main():
    ro.main()
