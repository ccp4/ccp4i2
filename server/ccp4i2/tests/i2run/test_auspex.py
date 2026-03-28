# Copyright (C) 2025 University of York
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


def test_auspex():
    args = ["AUSPEX"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    with i2run(args) as job:
        for plot in ["F", "FSigF", "I", "ISigI", "SigF", "SigI"]:
            assert (job / f"{plot}_plot.png").exists()
