# Copyright (C) 2026 University of York
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
import gemmi


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["SubtractNative"]
    args += ["--XYZIN", cif8xfm]
    args += ["--MAPIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        gemmi.read_ccp4_map(str(job / "MAPOUT.map"))
