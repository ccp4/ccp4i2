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
from shutil import which
import gemmi
from pytest import mark
from .utils import i2run


@mark.skipif(which("zanuda") is None and which("zanuda.exe") is None, reason="zanuda not installed")
def test_zanuda_8xfm(cif8xfm, mtz8xfm):
    """Test zanuda space group validation with 8xfm mmCIF input."""
    args = ["zanuda"]
    args += ["--XYZIN", f"fullPath={cif8xfm}"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]

    with i2run(args) as job:
        # Zanuda writes output as zanuda.pdb / zanuda.mtz (not XYZOUT)
        xyzout = job / "zanuda.pdb"
        assert xyzout.exists(), f"No zanuda.pdb: {list(job.iterdir())}"
        gemmi.read_pdb(str(xyzout))

        # Check split map coefficients
        for name in ["FPHIOUT", "DIFFPHIOUT"]:
            mtz_file = job / f"{name}.mtz"
            assert mtz_file.exists(), f"No {name}: {list(job.iterdir())}"
