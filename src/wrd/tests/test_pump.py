import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.pump import main


@pytest.mark.component
def test_pump_main():

    # August 19, 2021 Data
    # Stage 1 is default
    m = main()
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.33924  # kWh/m3
    # Stage 2
    m = main(Qin=1029, Pin=131.2 * pyunits.psi, stage_num=2)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.07753
    # Stage 3
    m = main(Qin=384, Pin=(112.6 - 41.9) * pyunits.psi, stage_num=3)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.40512
