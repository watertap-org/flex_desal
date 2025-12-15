import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.ro_stage import main

@pytest.mark.component
def test_ro_stage_main():
    m = main()