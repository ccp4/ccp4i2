"""Add to the GammaXe project (scenario_ep.py) the runs the Crank2 and SHELX
pages are illustrated from: the same xenon SAD data, solved from scratch,
the substructure found rather than given.

SHELX is licensed separately; tasks find it through the SHELXDIR preference
(or the environment variable of that name):

    env CCP4I2_HOME=/tmp/docs-home SHELXDIR=/path/to/shelx/bin \\
        ccp4-python ../docs/user/tools/scenario_ep_shelx.py

Adds: Crank2 (its full default pipeline) and SHELX (SHELXC/D/E, building,
refinement), then an unrun clone of each.
"""
import os
import sys
from pathlib import Path

from scenario_common import clone_last, i2run as _i2run, scratch_home
from scenario_ep import DEMO, PROJECT


def i2run(*args):
    _i2run(PROJECT, *args)


def main():
    if not os.environ.get("SHELXDIR"):
        sys.exit("Set SHELXDIR to the directory holding shelxc, shelxd, shelxe.")
    common = ["--F_SIGFanom", "fileOut=import_merged[-1].OBSOUT",
              "--SEQIN", f"fullPath={DEMO / 'gamma.asu.xml'}",
              "--ATOM_TYPE", "Xe", "--NUMBER_SUBSTRUCTURE", "2",
              "--WAVELENGTH", "1.54179", "--FPRIME", "-0.79", "--FDPRIME", "7.36"]
    i2run("crank2", *common)
    i2run("shelx", *common)
    for task in ("crank2", "shelx"):
        clone_last(PROJECT, task)
    # The SHELX run must be the SHELX route (it once ran crank2's).
    jobs = scratch_home() / "projects/gammaxe/CCP4_JOBS"
    logs = [p / "log.txt" for p in jobs.glob("job_*") if (p / "log.txt").exists()]
    assert any("Running shelxe" in log.read_text(errors="replace") for log in logs), \
        "no job ran SHELXE"

if __name__ == "__main__":
    main()
