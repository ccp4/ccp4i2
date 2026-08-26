import pytest
from .utils import i2run


@pytest.mark.skip(reason="NOT IMPLEMENTED: no test body. Would need DIALS-integrated serial .refl/.expt files; see the docstring.")
def test_ssx_reduce():
    """Test xia2.ssx_reduce serial data reduction.

    Requires DIALS-integrated serial crystallography reflection
    files (.refl) with matching experiment files (.expt).
    """
    pass
