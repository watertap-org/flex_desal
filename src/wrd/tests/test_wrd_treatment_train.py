import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main


@pytest.mark.unit
def test_no_wrd_file():
    with pytest.raises(
        ValueError, match="Input file must be provided to build WRD system."
    ):
        _ = main(num_pro_trains=2, file=None)


# Add test(s) in this file to simply check that the model builds and solves.
