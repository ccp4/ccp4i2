import pytest
from .utils import i2run


@pytest.mark.skip(reason="NOT IMPLEMENTED: no test body. Would need DIALS-integrated .refl/.expt files; see the docstring.")
def test_multiplex():
    """Test xia2.multiplex multi-crystal data processing.

    Requires multiple DIALS-integrated reflection files (.refl) with
    matching experiment files (.expt).
    """
    pass
