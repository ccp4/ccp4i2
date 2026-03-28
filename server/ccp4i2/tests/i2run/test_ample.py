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


def test_8xfm(mtz8xfm, seq8xfm):
    args = ["AMPLE"]
    args += ["--AMPLE_F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--AMPLE_SEQIN", seq8xfm]
    args += ["--AMPLE_EXISTING_MODELS", "False"]
    args += ["--AMPLE_NPROC", "8"]
    i2run(args)
