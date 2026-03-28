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
from .utils import demoData, i2run


def test_mdm2_two_models():
    """Test phaser.ensembler with two different MDM2 structures."""
    args = ["phaser_ensembler"]
    args += ["--XYZIN_LIST", demoData("mdm2", "4qo4.pdb")]
    args += ["--XYZIN_LIST", demoData("mdm2", "4hg7.pdb")]
    with i2run(args) as job:
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"
