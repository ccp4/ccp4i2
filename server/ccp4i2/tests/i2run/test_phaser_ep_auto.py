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
import gemmi
import pytest
from .utils import demoData, i2run


# Run phaser tests first to avoid RDKit pickle contamination
@pytest.mark.order("first")
def test_gamma_xe():
    """Test phaser_EP_AUTO with gamma Xe anomalous data."""
    args = ["phaser_EP_AUTO"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_HA", demoData("gamma", "heavy_atoms.pdb")]
    args += ["--ASUFILE", demoData("gamma", "gamma.asu.xml")]
    args += ["--COMP_BY", "ASU"]
    args += ["--PARTIALMODELORMAP", "NONE"]
    args += ["--WAVELENGTH", "1.542"]
    args += ["--LLGC_CYCLES", "20"]
    args += ["--ELEMENTS", "Xe"]

    with i2run(args) as job:
        # Check output PDB files exist
        pdb_files = list(job.glob("PHASER*.pdb"))
        assert len(pdb_files) >= 1, f"No Phaser PDB output: {list(job.iterdir())}"

        # Check output MTZ files
        mtz_files = list(job.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output: {list(job.iterdir())}"
