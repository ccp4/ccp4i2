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


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["sheetbend"]
    args += ["--XYZIN", cif8xfm]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
        xml = ET.parse(job / "program.xml")
        rwork = float(xml.find(".//Final/Rwork").text)
        rfree = float(xml.find(".//Final/Rfree").text)
        assert rwork < 0.23
        assert rfree < 0.25
