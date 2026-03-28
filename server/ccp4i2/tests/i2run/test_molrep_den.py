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
from .urls import redo_mtz
from .utils import demoData, download, i2run


def test_gamma_density():
    """Test molrep density-mode search with PDB-REDO map coefficients."""
    with download(redo_mtz("1gyv")) as redo:
        args = ["molrep_den"]
        args += ["--F_PHI_MAP", f"fullPath={redo}", "columnLabels=/*/*/[FWT,PHWT]"]
        args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
        with i2run(args) as job:
            xyzout = job / "XYZOUT.pdb"
            assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
            st = gemmi.read_structure(str(xyzout))
            assert len(st[0]) > 0
