import pytest
from .utils import i2run


@pytest.mark.skip(reason="Requires CrystFEL serial crystallography input data not in demo_data")
def test_serial_import():
    """Test import_serial_pipe serial crystallography data import.

    This pipeline runs import_serial (CrystFEL→MTZ conversion) then
    aimless_pipe for analysis.  Needs CrystFEL .hkl half-dataset files.
    """
    pass
