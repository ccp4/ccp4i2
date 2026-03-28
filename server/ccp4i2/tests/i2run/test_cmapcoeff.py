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
from .utils import demoData, i2run


def test_gamma_anomalous():
    """Test anomalous map coefficient calculation with gamma Xe data."""
    mtz = demoData("gamma", "merged_intensities_Xe.mtz")
    args = ["cmapcoeff"]
    args += ["--MAPTYPE", "anom"]
    args += ["--F_SIGF1", mtz]
    args += ["--ABCD1", demoData("gamma", "initial_phases.mtz")]
    with i2run(args) as job:
        fphiout = job / "FPHIOUT.mtz"
        assert fphiout.exists(), f"No FPHIOUT: {list(job.iterdir())}"
        gemmi.read_mtz_file(str(fphiout))
