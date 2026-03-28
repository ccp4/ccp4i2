# Copyright (C) 2025-2026 Newcastle University
# Copyright (C) 2025-2026 University of York
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
from .utils import i2run


def test_lorestr_1h1s(pdb1h1s, mtz1h1s):
    """Test lorestr low-resolution refinement with 1h1s data (no auto homologue fetch)."""
    args = ["lorestr_i2"]
    args += ["--XYZIN", f"fullPath={pdb1h1s}"]
    args += ["--F_SIGF", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FREE]"]
    args += ["--AUTO", "none"]
    args += ["--CPU", "1"]

    with i2run(args) as job:
        # Check refined model output
        xyzout = job / "XYZOUT.pdb"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
        gemmi.read_pdb(str(xyzout))

        # Check map coefficients
        for name in ["FPHIOUT", "DIFFPHIOUT"]:
            mtz_file = job / f"{name}.mtz"
            assert mtz_file.exists(), f"No {name}: {list(job.iterdir())}"


def test_lorestr_8xfm_mmcif(cif8xfm, mtz8xfm):
    """Test lorestr with mmCIF input containing long ligand residue codes (8xfm)."""
    args = ["lorestr_i2"]
    args += ["--XYZIN", f"fullPath={cif8xfm}"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    args += ["--AUTO", "none"]
    args += ["--CPU", "1"]

    with i2run(args) as job:
        # Check refined model output
        xyzout = job / "XYZOUT.pdb"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"

        # Check map coefficients
        for name in ["FPHIOUT", "DIFFPHIOUT"]:
            mtz_file = job / f"{name}.mtz"
            assert mtz_file.exists(), f"No {name}: {list(job.iterdir())}"
