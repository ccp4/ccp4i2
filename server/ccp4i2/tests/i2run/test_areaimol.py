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
import gemmi

from .urls import rcsb_pdb
from .utils import download, i2run


def test_5a3h():
    with download(rcsb_pdb("5a3h")) as pdb_5a3h:
        args = ["areaimol"]
        args += ["--XYZIN", pdb_5a3h]
        args += ["--OUTPUT_MODE", "ATOM"]
        with i2run(args) as job:
            gemmi.read_structure(str(job / "XYZOUT.pdb"))


def test_5a3h_imol():
    with download(rcsb_pdb("5a3h")) as pdb_5a3h:
        args = ["areaimol"]
        args += ["--XYZIN", pdb_5a3h]
        args += ["--DIFFMODE", "IMOL"]
        args += ["--OUTPUT_MODE", "ATOM"]
        with i2run(args) as job:
            gemmi.read_structure(str(job / "XYZOUT.pdb"))


def test_compare_5a3h_and_8a3h():
    with download(rcsb_pdb("5a3h")) as pdb_5a3h, download(rcsb_pdb("8a3h")) as pdb_8a3h:
        args = ["areaimol"]
        args += ["--XYZIN", pdb_5a3h]
        args += ["--XYZIN2", pdb_8a3h]
        args += ["--DIFFMODE", "COMPARE"]
        args += ["--OUTPUT_MODE", "ATOM"]
        with i2run(args) as job:
            gemmi.read_structure(str(job / "XYZOUT.pdb"))


# def test_8xfm(cif8xfm):
#     args = ["areaimol"]
#     args += ["--XYZIN", cif8xfm]
#     args += ["--OUTPUT_MODE", "ATOM"]
#     with i2run(args) as job:
#         assert hasLongLigandName(job / "XYZOUT.cif")
