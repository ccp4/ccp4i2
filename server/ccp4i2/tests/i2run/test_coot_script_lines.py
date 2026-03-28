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
import pytest
from .utils import demoData, i2run


@pytest.mark.skip(reason="Requires coot --no-graphics binary; may not be available in CI")
def test_rename_chain():
    """Test scripted model building with coot."""
    args = ["coot_script_lines"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--SCRIPT", "rename_chain(0, 'A', 'B')"]
    with i2run(args) as job:
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"
