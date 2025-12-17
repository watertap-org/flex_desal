import pytest
from wrd.components.uf_system import main


@pytest.mark.component
def test_uf_system_with_costing():
    m = main(add_costing=True)
