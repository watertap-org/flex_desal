import pytest

from pyomo.environ import value, assert_optimal_termination, units as pyunits

from watertap.core.solvers import get_solver

import srp.components.pump as pump


@pytest.mark.component
def test_pump_main():
    pump.main()


@pytest.mark.component
def test_pump_operation():
    m = pump.build_system()
    pump.set_system_scaling(m)
    pump.set_system_op_conditions(m)
    pump.set_pump_op_conditions(m.fs.pump, pressure=400, efficiency=0.85)
    pump.init_system(m)
    solver = get_solver()
    results = solver.solve(m.fs)
    assert_optimal_termination(results)

    assert value(m.fs.pump.unit.work_mechanical[0]) == pytest.approx(
        2236626.4, rel=1e-3
    )
