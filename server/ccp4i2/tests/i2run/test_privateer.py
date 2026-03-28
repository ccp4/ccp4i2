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


def test_glyco_4iid():
    """Test privateer glycan validation with 4iid."""
    args = ["privateer"]
    args += ["--XYZIN", demoData("glyco", "4iid.pdb")]
    args += ["--F_SIGF", f"fullPath={demoData('glyco', '4iid.mtz')}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--SINGLETHREADED", "True"]

    with i2run(args) as job:
        assert (job / "program.xml").exists(), "No program.xml output"
