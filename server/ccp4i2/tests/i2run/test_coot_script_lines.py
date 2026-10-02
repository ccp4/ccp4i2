"""Scripted Coot, run headless with CCP4's Coot 1.

It was skipped unconditionally, and broken twice over: the wrapper launched
"coot", which a CCP4 9 install does not have (coot-1), and Coot 1's
scripting functions are in modules, not the script's globals.
"""
import shutil

import gemmi
import pytest

from .utils import demoData, i2run

pytestmark = pytest.mark.skipif(shutil.which("coot-1") is None, reason="needs CCP4's coot-1")


def test_rename_chain():
    args = ["coot_script_lines"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--SCRIPT", "change_chain_id(MolHandle_1, 'A', 'B', 0, 0, 0)\n"
                         "write_pdb_file(MolHandle_1, os.path.join(dropDir, 'output.pdb'))\n"]
    with i2run(args) as job:
        out = [p for p in job.rglob("*.pdb")]
        assert out, f"No PDB output: {list(job.iterdir())}"
        chains = {c.name for p in out for c in gemmi.read_structure(str(p))[0]}
        assert "B" in chains, chains
        # Named by the recipe, not only the file the script wrote.
        from ccp4i2.db import models
        record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
        names = [f.annotation for f in models.File.objects.filter(
            job=record, job_param_name__startswith="XYZOUT")]
        assert names == ["Scripted Coot, own script: output.pdb"], names
