import pytest
from pyomo.environ import units as pyunits
from wrd.components.pump import main


@pytest.mark.component
def test_pump_main():
    m = main()
    m = main(Qin=1029, Pin=131.2 * pyunits.psi, stage_num=2)
    m = main(Qin=384, Pin=112.6 * pyunits.psi, stage_num=3)
