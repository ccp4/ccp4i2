import pytest
from .utils import i2run


@pytest.mark.skip(reason="Requires DIALS integrated .refl/.expt files not in demo_data")
def test_multiplex():
    """Test xia2.multiplex multi-crystal data processing.

    Requires multiple DIALS-integrated reflection files (.refl) with
    matching experiment files (.expt).
    """
    pass
