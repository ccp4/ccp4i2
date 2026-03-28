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
import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


def test_dimple():
    args = ["i2Dimple"]  # Plugin is named i2Dimple in registry
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_native.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["FPHIOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
