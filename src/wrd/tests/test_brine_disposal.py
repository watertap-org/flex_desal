import pytest

from wrd.components.brine_disposal import main


@pytest.mark.component
def test_brine_disposal():
    _ = main()
