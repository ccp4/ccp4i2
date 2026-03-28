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
import json
from pytest import approx
import gemmi
from .utils import demoData, i2run


def test_mrparse():
    args = ["mrparse"]
    args += ["--SEQIN", demoData("gamma", "gamma.pir")]
    args += ["--DATABASE", "PDB"]
    args += ["--USEAPI", "False"]
    with i2run(args) as job:
        outputs = list(job.glob("*.pdb"))
        assert len(outputs) > 0
        for output in outputs:
            gemmi.read_structure(str(output))
        with (job / "mrparse_0" / "homologs.json").open() as f:
            homologues = json.load(f)
            assert max(h["seq_ident"] for h in homologues) == approx(1.0)
