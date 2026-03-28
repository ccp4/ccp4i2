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


@pytest.mark.skip(reason="Requires DIALS integrated .refl/.expt files not in demo_data")
def test_multiplex():
    """Test xia2.multiplex multi-crystal data processing.

    Requires multiple DIALS-integrated reflection files (.refl) with
    matching experiment files (.expt).
    """
    pass
