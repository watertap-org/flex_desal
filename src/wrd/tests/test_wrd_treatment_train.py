import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main


@pytest.mark.component  # I don't have a good sense of unit vs. component
def test_no_wrd_file():
    with pytest.raises(
        ValueError, match="Input file must be provided to build WRD system."
    ):
        _ = main(num_pro_trains=2, file=None)


# @pytest.mark.parametrize("num_pro_trains", [1,2,3,4])
# @pytest.mark.component
# def test_wrd_treatment_train(num_pro_trains):
#     file = None
#     Qin = (num_pro_trains * 2600) * pyunits.m**3 / pyunits.day
#     m = main(num_pro_trains=num_pro_trains, uf_split_fraction=[0.4, 0.4, 0.2], file=file)
# # Add test(s) in this file to simply check that the model builds and solves.
