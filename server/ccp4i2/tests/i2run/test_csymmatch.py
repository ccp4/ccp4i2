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
from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm):
    args = ["csymmatch"]
    args += ["--XYZIN_QUERY", cif8xfm]
    args += ["--XYZIN_TARGET", cif8xfm]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
        xml = ET.parse(job / "program.xml")
        change = xml.find(".//ChangeOfOrigin").text.replace(" ", "")
        assert change == "uvw=(0,0,0)"
