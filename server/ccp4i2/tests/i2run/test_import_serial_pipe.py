import pytest
from .utils import i2run


@pytest.mark.skip(reason="NOT IMPLEMENTED: no test body. Would need CrystFEL .hkl half-datasets; see the docstring.")
def test_serial_import():
    """Test import_serial_pipe serial crystallography data import.

    This pipeline runs import_serial (CrystFEL→MTZ conversion) then
    aimless_pipe for analysis.  Needs CrystFEL .hkl half-dataset files.
    """
    pass
