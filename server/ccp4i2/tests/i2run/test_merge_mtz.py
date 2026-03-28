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
import gemmi
from .utils import demoData, i2run


@pytest.mark.skip(reason="CLI list parsing creates CString items instead of CMergeMiniMtzItem — wrapper expects .fileName attribute")
def test_gamma_merge_two_mtz():
    """Test merging two mini-MTZ files from gamma data."""
    mtz1 = demoData("gamma", "merged_intensities_native.mtz")
    mtz2 = demoData("gamma", "merged_intensities_Xe.mtz")
    args = ["mergeMtz"]
    args += ["--MINIMTZINLIST", f"filename={mtz1}"]
    args += ["--MINIMTZINLIST", f"filename={mtz2}"]
    with i2run(args) as job:
        hklout = job / "HKLOUT.mtz"
        assert hklout.exists(), f"No HKLOUT: {list(job.iterdir())}"
        gemmi.read_mtz_file(str(hklout))
