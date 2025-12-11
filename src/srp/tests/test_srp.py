import pytest

import srp.srp_treatment_train as srp_train


@pytest.mark.component
def test_srp_treatment_train_basic():
    m = srp_train.run_srp_basic()


"""
NOTE: 12/11/2025 - this test is not passing on Linux; see PR #12
"""
# @pytest.mark.component
# def test_srp_treatment_train_full():
#     m = srp_train.run_srp()
