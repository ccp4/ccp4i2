import pytest
from .utils import i2run


@pytest.mark.skip(reason="Requires XIA2 output directory structure not in demo_data")
def test_import_xia2():
    """Test AlternativeImportXIA2 harvesting of xia2 results.

    This task scans a xia2 output directory for integrated/merged
    MTZ files, ispyb.xml, and log files.  Needs a real xia2 run directory.
    """
    pass
