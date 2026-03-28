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


@pytest.mark.skip(reason="Requires external PDB-REDO web service; not suitable for offline CI")
def test_pdb_redo():
    """Test PDB-REDO web services task.

    Submits a structure to the PDB-REDO server for re-refinement.
    Requires network access and may take several minutes.
    """
    pass
