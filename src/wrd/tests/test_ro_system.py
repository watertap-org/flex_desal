import pytest
from wrd.components.ro_system import main


@pytest.mark.component
def test_ro_system_with_costing():
    _ = main(add_costing=True)


@pytest.mark.component
def test_ro_system_without_costing():
    _ = main(add_costing=False)


# Not sure what other tests would be useful here
