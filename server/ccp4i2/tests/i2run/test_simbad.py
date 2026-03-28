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
from .utils import demoData, i2run


@pytest.mark.slow
def test_gamma_lattice():
    """Test SIMBAD lattice search with gamma native data.

    Marked slow because SIMBAD lattice search queries a database.
    """
    args = ["SIMBAD"]
    args += ["--F_SIGF", demoData("gamma", "gamma_native.mtz")]
    args += ["--SIMBAD_SEARCH_LEVEL", "Lattice"]
    args += ["--SIMBAD_NPROC", "1"]
    with i2run(args) as job:
        assert (job / "diagnostic.xml").exists()
