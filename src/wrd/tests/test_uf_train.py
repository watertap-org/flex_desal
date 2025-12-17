import pytest
from wrd.components.uf_train import main


@pytest.mark.component
def test_uf_train_with_costing():
    m = main(add_costing=True)
