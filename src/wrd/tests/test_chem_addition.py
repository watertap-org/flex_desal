import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_equivalent
import wrd.components.chemical_addition as ca


@pytest.mark.unit
def test_chem_addition_missing_data():
    msg = "Must specify a chemical for addition."
    with pytest.raises(ValueError, match=msg):
        _ = ca.build_system(
            chemical_name=None
        )  # Why is this build system and not main?


@pytest.mark.component
def test_chem_addition():
    _ = ca.main()


@pytest.mark.component
def test_chem_addition_setting_dose():
    # Test that component can accept dose and chem cost directly instead of reading from yaml

    expected_chem_flow = 150.1  # gal/day
    # = 10000 gal/min * (60*24 min/day) * 5 mg/L * (1 L / 0.264172 gal) * (1 kg / 1e6 mg) * (1 m3 / 1200 kg) * (264.172 gal/m3) / (0.4 gal_chem/gal_sol)
    expected_chem_cost = 300.2  # $/day
    # = 150.1 gal/day * 2 $/gal = 300.2 $/day
    m = ca.main(chemical_name="sodium_bisulfite", Qin=10000, dose=5, chem_cost=2)
    # Note that ratio in solution and solution density are still read from yaml!
    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.chem_addition.unit.chemical_soln_flow_vol,
                    to_units=pyunits.gallon / pyunits.day,
                )
            ),
            rel=1e-3,
        )
        == expected_chem_flow
    )
    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.costing.aggregate_flow_costs["sodium_bisulfite"],
                    to_units=pyunits.USD_2021 / pyunits.day,
                )
            ),
            rel=1e-3,
        )
        == expected_chem_cost
    )


# Could add test for passing different costing units. Right now, only USD_2021 / gallon is supported.
