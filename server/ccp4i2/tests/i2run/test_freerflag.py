# Copyright (C) 2025 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
from pytest import approx
import gemmi
from .utils import demoData, i2run


def test_freerflag():
    args = ["freerflag"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--COMPLETE", "False"]
    args += ["--FRAC", "0.2"]
    with i2run(args) as job:
        mtz = gemmi.read_mtz_file(str(job / "FREEROUT.mtz"))
        column = mtz.rfree_column()
        fraction = list(column).count(0) / len(column)
        assert fraction == approx(0.2, abs=0.01)
