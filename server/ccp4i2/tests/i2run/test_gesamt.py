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
import xml.etree.ElementTree as ET
import gemmi
from .utils import i2run


def test_self_superposition(cif8xfm):
    """Test gesamt pairwise superposition (identity)."""
    args = ["gesamt"]
    args += ["--XYZIN_QUERY", cif8xfm]
    args += ["--XYZIN_TARGET", cif8xfm]

    with i2run(args) as job:
        gemmi.read_structure(str(job / "XYZOUT.cif"))

        # Self-superposition should give RMSD ~ 0
        tree = ET.parse(job / "program.xml")
        rms = tree.find(".//Transformation")
        if rms is not None:
            assert float(rms.get("rms", "1")) < 0.1
