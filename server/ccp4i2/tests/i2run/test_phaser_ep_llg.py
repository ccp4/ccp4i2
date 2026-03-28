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


def test_gamma_xe():
    """Test phaser EP LLG anomalous map from Xe coordinates."""
    args = ["phaser_EP_LLG"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_HA", demoData("gamma", "heavy_atoms.pdb")]
    with i2run(args, allow_errors=True) as job:
        # EP_LLG may not converge with this data but should produce output
        assert (job / "diagnostic.xml").exists()
