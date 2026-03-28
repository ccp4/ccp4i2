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
from .utils import demoData, i2run


def test_gamma_phase_comparison():
    """Test phase comparison between two sets of phases."""
    args = ["cphasematch"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_native.mtz")]
    args += ["--ABCD1", demoData("gamma", "initial_phases.mtz")]
    args += ["--ABCD2", demoData("gamma", "initial_phases.mtz")]
    with i2run(args) as job:
        abcdout = job / "ABCDOUT.mtz"
        assert abcdout.exists(), f"No ABCDOUT: {list(job.iterdir())}"
