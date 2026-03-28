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


@pytest.mark.skip(reason="Requires CrystFEL serial crystallography input data not in demo_data")
def test_serial_import():
    """Test import_serial_pipe serial crystallography data import.

    This pipeline runs import_serial (CrystFEL→MTZ conversion) then
    aimless_pipe for analysis.  Needs CrystFEL .hkl half-dataset files.
    """
    pass
