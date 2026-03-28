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
from .utils import demoData, i2run


def test_gamma_model():
    """Test density calculation from gamma model."""
    args = ["density_calculator"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--D_MIN", "2.0"]
    with i2run(args) as job:
        mapout = job / "MAPOUT.map"
        assert mapout.exists(), f"No MAPOUT: {list(job.iterdir())}"
        ccp4map = gemmi.read_ccp4_map(str(mapout))
        assert ccp4map.grid.nu > 0
