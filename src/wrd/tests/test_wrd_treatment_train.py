import pytest
from pyomo.environ import value, units as pyunits
from wrd.wrd_treatment_train import main


@pytest.mark.unit
def test_no_wrd_file():
    with pytest.raises(
        ValueError, match="Input file must be provided to build WRD system."
    ):
        _ = main(num_pro_trains=2, file=None)


@pytest.mark.parametrize("num_pro_trains", [1, 2, 3, 4])
@pytest.mark.component
def test_wrd_treatment_train_8_19_21(num_pro_trains):
    file = "wrd_inputs_8_19_21.yaml"
    m = main(num_pro_trains=num_pro_trains, file=file)
    # TODO: Add tests against facility data


@pytest.mark.parametrize("num_pro_trains", [1, 2, 3, 4])
@pytest.mark.component
def test_wrd_treatment_train_3_13_21(num_pro_trains):
    file = "wrd_inputs_3_13_21.yaml"
    m = main(num_pro_trains=num_pro_trains, file=file)
    # TODO: Add tests against facility data
