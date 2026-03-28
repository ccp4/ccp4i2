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
from .utils import i2run


def test_8xfm(cif8xfm, mtz8xfm):
    """Test edstats with 8xfm structure and map coefficients."""
    args = ["edstats"]
    args += ["--XYZIN", cif8xfm]
    args += ["--FPHIIN1", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    args += ["--FPHIIN2", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[DELFWT,PHDELWT]"]
    args += ["--RES_LOW", "50.0"]
    args += ["--RES_HIGH", "1.3"]

    with i2run(args) as job:
        assert (job / "program.xml").exists(), "No program.xml output"
