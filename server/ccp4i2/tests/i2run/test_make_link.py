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
from gemmi import cif, read_pdb
from .urls import pdbe_pdb
from .utils import download, i2run


def test_6ndn():
    with download(pdbe_pdb("6ndn")) as pdb:
        args = ["MakeLink"]
        args += ["--RES_NAME_1_TLC", "LYS"]
        args += ["--RES_NAME_2_TLC", "PLP"]
        args += ["--ATOM_NAME_1_TLC", "NZ"]
        args += ["--ATOM_NAME_2_TLC", "C4A"]
        args += ["--ATOM_NAME_1", "NZ"]
        args += ["--ATOM_NAME_2", "C4A"]
        args += ["--TOGGLE_DELETE_2", "True"]
        args += ["--DELETE_2", "O4A"]
        args += ["--BOND_ORDER", "DOUBLE"]
        args += ["--TOGGLE_LINK", "True"]
        args += ["--XYZIN", pdb]
        # allow_errors=True: Allow diagnostic warnings for empty optional strings
        # (CHARGE_1_LIST, DELETE_1_LIST, etc.) - the script handles empty values correctly
        with i2run(args, allow_errors=True) as job:
            doc = cif.read(str(job / "LYS-PLP_link.cif"))
            for name in ("mod_LYSm1", "mod_PLPm1", "link_LYS-PLP"):
                assert name in doc
            read_pdb(str(job / "ModelWithLinks.pdb"))
