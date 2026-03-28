# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
import pytest
from .utils import i2run


@pytest.mark.skip(reason="Requires XIA2 output directory structure not in demo_data")
def test_import_xia2():
    """Test AlternativeImportXIA2 harvesting of xia2 results.

    This task scans a xia2 output directory for integrated/merged
    MTZ files, ispyb.xml, and log files.  Needs a real xia2 run directory.
    """
    pass
