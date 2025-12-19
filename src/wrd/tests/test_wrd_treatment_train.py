import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main


@pytest.mark.parametrize("num_pro_trains", [4, 3, 2])
@pytest.mark.component
def test_wrd_treatment_train(num_pro_trains):
    m = main(num_pro_trains=num_pro_trains)
