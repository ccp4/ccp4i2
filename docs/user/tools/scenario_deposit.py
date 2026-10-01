"""Add to the MDM2 project the runs three more pages are illustrated from:

- Pdbset, a scripted edit: the refined model cut back to poly-alanine
  (EXCLUDE SIDE), the kind of search model molecular replacement uses when
  the side chains are not to be trusted, with its chain renamed.
- SubtractNative: the refined model's calculated density, half of it,
  taken from the 2mFo-DFc map, leaving density the model does not explain.
- Preparing the refined model for deposition. Its inputs are found the way
  the task's interface finds them, by tracing the model back through the
  refinement and scaling jobs (trace_xyzin_lineage). It is NOT sent to the
  wwPDB validation server: that uploads the structure to an outside service.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_deposit.py

Run scenario_refine.py, scenario_tools.py and scenario_refine_tools.py first.
Each task's last job gets an unrun clone.
"""
import json
import os
import subprocess
import sys

from scenario_common import DJANGO, clone_last, i2run as _i2run
from scenario_refine import PROJECT


def i2run(*args):
    _i2run(PROJECT, *args)


def deposition_inputs() -> dict:
    """What the deposition task's interface fills in from the refined model."""
    out = subprocess.run([sys.executable, "manage.py", "shell", "-c", (
        "import json\n"
        "from ccp4i2.db.models import File\n"
        "from ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.adding_stats_to_mmcif_i2 "
        "import adding_stats_to_mmcif_i2 as task\n"
        # The refinement's mmCIF model, which carries its statistics; the
        # PDB-format one cannot be deposited and the task refuses it.
        f"f = File.objects.filter(job__project__name='{PROJECT}', job__task_name='prosmart_refmac', "
        "job_param_name='CIFFILE', job__parent__isnull=True).order_by('job__id').first()\n"
        "r = task.trace_xyzin_lineage(task.__new__(task), str(f.uuid))\n"
        "r['xyzin'] = str(f.uuid)\n"
        "print('TRACE=' + json.dumps(r, default=str))\n")],
        check=True, env={**os.environ, **DJANGO}, capture_output=True, text=True).stdout
    trace = json.loads(next(line[6:] for line in out.splitlines() if line.startswith("TRACE=")))
    assert trace["success"], trace
    return trace


def main():
    refined = "fileOut=prosmart_refmac[-1].XYZOUT"

    i2run("pdbset_ui", "--XYZIN", refined,
          "--EXTRA_PDBSET_KEYWORDS", "EXCLUDE SIDE\nCHAIN M")

    i2run("SubtractNative", "--MAPIN", "fileOut=prosmart_refmac[-1].FPHIOUT",
          "--XYZIN", refined, "--FRACTION", "0.5")

    trace = deposition_inputs()
    args = ["--XYZIN", f"dbFileId={trace['xyzin']}",
            "--ASUCONTENT", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE",
            "--SENDTOVALIDATIONSERVER", "False"]
    for param, file_id in trace["files"].items():
        args += [f"--{param}", f"dbFileId={file_id}"]
    for param, path in trace["paths"].items():
        args += [f"--{param}", f"fullPath={path}"]
    for param, value in trace["params"].items():
        args += [f"--{param}", str(value)]
    i2run("adding_stats_to_mmcif_i2", *args)

    for task in ("pdbset_ui", "SubtractNative", "adding_stats_to_mmcif_i2"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
