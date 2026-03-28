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
from sys import platform
from gemmi import read_mtz_file, read_pdb
from pytest import fixture, mark
from .urls import redo_mtz
from .utils import download, i2run


pytestmark = mark.skipif(platform == "win32", reason="Not supported on Windows")


@fixture(scope="module", name="mtz")
def mtz_fixture():
    with download(redo_mtz("2ccf")) as path:
        yield path


@mark.skip(reason="Skipping temporarily as no shelx in build")
def test_arcimboldo(mtz):
    args = ["arcimboldo"]
    args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--ARCIMBOLDO_OPTIONS", "LITE"]
    args += ["--N_COMPONENTS", "1"]
    args += ["--MOLECULAR_WEIGHT", "6000"]
    args += ["--N_FRAGMENTS", "1"]
    args += ["--HELIX_LENGTH", "14"]
    with i2run(args) as job:
        read_mtz_file(str(job / "anis.mtz"))
        read_pdb(str(job / "anis.patterson.pdb"))
        log = (job / "log.txt").read_text()
        assert "TNCS was found. Please add another fragment" in log
