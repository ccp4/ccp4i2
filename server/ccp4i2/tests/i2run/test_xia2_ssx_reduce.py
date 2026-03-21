import pytest
from .utils import i2run


@pytest.mark.skip(reason="Requires DIALS SSX integrated .refl/.expt files not in demo_data")
def test_ssx_reduce():
    """Test xia2.ssx_reduce serial data reduction.

    Requires DIALS-integrated serial crystallography reflection
    files (.refl) with matching experiment files (.expt).
    """
    pass
