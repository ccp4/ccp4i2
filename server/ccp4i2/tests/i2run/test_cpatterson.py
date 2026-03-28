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


def test_gamma_native():
    """Test Patterson map calculation from gamma native data."""
    args = ["cpatterson"]
    args += ["--F_SIGF", demoData("gamma", "gamma_native.mtz")]
    with i2run(args) as job:
        assert (job / "MAPOUT.map").exists(), f"No MAPOUT: {list(job.iterdir())}"
