import pytest
from .utils import i2run


@pytest.mark.skip(reason="NOT IMPLEMENTED: no test body. Would need a real xia2 output directory; see the docstring.")
def test_import_xia2():
    """Test AlternativeImportXIA2 harvesting of xia2 results.

    This task scans a xia2 output directory for integrated/merged
    MTZ files, ispyb.xml, and log files.  Needs a real xia2 run directory.
    """
    pass
