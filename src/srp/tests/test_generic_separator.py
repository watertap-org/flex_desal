import pytest

from pyomo.environ import value, assert_optimal_termination, units as pyunits

from watertap.core.solvers import get_solver

import srp.components.generic_separator as separator


@pytest.mark.component
def test_generic_separator_main():
    separator.main()


@pytest.mark.component
def test_generic_separator_outlets():
    m = separator.build_system(
        Qin=10000, Cin=1200, outlet_list=["cats", "dogs", "fish"]
    )
    assert len(m.fs.separator.unit.config.outlet_list) == 3
    for i, outlet in enumerate(m.fs.separator.unit.config.outlet_list, 1):
        prod_unit = m.fs.find_component(f"product{i}")
        assert prod_unit is not None
    separator.set_system_scaling(m)
    separator.set_system_op_conditions(m)
    separator.init_system(m)
    solver = get_solver()
    results = solver.solve(m.fs)
    assert_optimal_termination(results)

    for outlet in m.fs.separator.unit.config.outlet_list:
        sep_out = m.fs.separator.find_component(f"{outlet}")
        assert value(
            pyunits.convert(
                sep_out.properties[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallon / pyunits.minute,
            )
        ) == pytest.approx(
            value(m.Qin) / len(m.fs.separator.unit.config.outlet_list), rel=1e-3
        )
        assert value(
            pyunits.convert(
                sep_out.properties[0].conc_mass_phase_comp["Liq", "TDS"],
                to_units=pyunits.mg / pyunits.L,
            )
        ) == pytest.approx(value(m.Cin), rel=1e-3)
