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


@pytest.mark.skip(reason="KeywordExtractor bug: get_merged_metadata is None for ctruncate container")
def test_gamma_intensities():
    """Test ctruncate with gamma anomalous intensities."""
    args = ["ctruncate"]
    args += ["--OBSIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--SEQIN", demoData("gamma", "gamma.pir")]

    with i2run(args) as job:
        # Check for output MTZ
        mtz_files = list(job.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output from ctruncate: {list(job.iterdir())}"
