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


def test_reindex_to_coords():
    """Test pointless reindexing gamma data to match coordinate reference."""
    args = ["pointless_reindexToMatch"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_REF", demoData("gamma", "gamma_model.pdb")]
    args += ["--REFERENCE", "XYZIN_REF"]

    with i2run(args) as job:
        mtz_files = list(job.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output from pointless: {list(job.iterdir())}"
