import pytest

from pyomo.environ import value, assert_optimal_termination, units as pyunits

from watertap.core.solvers import get_solver

import srp.components.generic_mixer as mixer


@pytest.mark.component
def test_generic_mixer_main():
    mixer.main()


@pytest.mark.component
def test_generic_mixer_inlets():
    m = mixer.build_system(inlet_list=["cats", "dogs", "fish"])
    assert len(m.inlet_list) == 3
    for inlet in m.inlet_list:
        feed_unit = m.fs.find_component(f"{inlet}")
        assert feed_unit is not None
    mixer.set_system_scaling(m)
    mixer.set_system_op_conditions(m)
    mixer.init_system(m)
    solver = get_solver()
    results = solver.solve(m.fs)
    assert_optimal_termination(results)

    for inlet in m.inlet_list:
        mix_in = m.fs.find_component(f"{inlet}")
        assert value(
            pyunits.convert(
                mix_in.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallon / pyunits.minute,
            )
        ) == pytest.approx(value(m.Qin) / len(m.inlet_list), rel=1e-3)
        assert value(
            pyunits.convert(
                mix_in.properties[0].conc_mass_phase_comp["Liq", "TDS"],
                to_units=pyunits.mg / pyunits.L,
            )
        ) == pytest.approx(value(m.Cin), rel=1e-3)
