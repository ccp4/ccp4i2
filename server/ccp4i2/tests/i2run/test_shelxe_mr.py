# Copyright (C) 2025-2026 Newcastle University
# Copyright (C) 2025-2026 University of York
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
from shutil import which
from gemmi import read_mtz_file, read_pdb
from pytest import mark
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


@mark.skipif(which("shelxe") is None and which("shelxe.exe") is None, reason="SHELXE not installed")
def test_gamma():
    args = ["shelxeMR"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--NTCYCLES", "2"]
    args += ["--NMCYCLES", "2"]
    with i2run(args) as job:
        read_pdb(str(job / "shelxrun.pdb"))
        read_mtz_file(str(job / "FPHIOUT.mtz"))
        xml = ET.parse(job / "program.xml")
        assert float(xml.find(".//BestCC").text) > 43
