import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.pump import main


@pytest.mark.skip  # Pending new pump model updates
def test_pump_main():
    # August 19, 2021 Data
    # Stage 1 is default
    m = main()
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.33924  # kWh/m3
    # Stage 2
    m = main(Qin=1029, Pin=131.2 * pyunits.psi, stage_num=2)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.07753
    # Stage 3
    m = main(Qin=384, Pin=112.6 * pyunits.psi, stage_num=3)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.22392


# TODO: Update these tests after the pump model updates have been implemented.
# Add uf pump test
# Removing numbers from the yaml that are only used for testing and hard coding them into the test file.
# Add tests for calculating pump speed and efficiency, especially for out of bounds cases.
