import pytest

import srp.srp_treatment_train as srp_train


@pytest.mark.component
def test_srp_treatment_train():
    m = srp_train.run_srp()
