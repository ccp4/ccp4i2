import pytest
from .utils import demoData, i2run


@pytest.mark.skip(reason="matthews has no Python wrapper script (def.xml only); needs investigation")
def test_gamma():
    """Test matthews coefficient calculation with gamma data."""
    pass
