"""Build the project the import tasks' page is illustrated from: each kind of
file brought into a project on its own.

Most come from the gamma demo (demo_data/gamma): a model, a sequence, AU
contents, unmerged data. A dictionary comes from the baz2b demo. The rest
come from one *full* MTZ file, the situation these import tasks exist for: a
file from another program that holds several kinds of data at once. It is the
complete output of a refinement in the Gamma project (scenario_maps.py), so
it has amplitudes, a free R set, map coefficients and phases together, and
each import extracts one kind. A map file is calculated from its map
coefficients with gemmi, as another program would have written it.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_import.py

Run scenario_maps.py first. Makes the project Imports; each task's last job
gets an unrun clone.
"""
import shutil
from pathlib import Path

import gemmi

from scenario_common import clone_last, i2run as _i2run, inputs_dir, scratch_home

PROJECT = "Imports"
DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data"
GAMMA = DEMO / "gamma"


def i2run(*args):
    _i2run(PROJECT, *args)


def full_mtz() -> Path:
    """The Gamma project's refinement output, copied under a name a user
    might have: refined_gamma.mtz."""
    found = sorted((scratch_home() / "projects/gamma/CCP4_JOBS").glob("job_*/COMPLETE_MTZ.mtz"),
                   key=lambda p: p.stat().st_mtime)
    if not found:
        raise SystemExit("No refinement in the Gamma project: run scenario_maps.py first.")
    target = inputs_dir() / "refined_gamma.mtz"
    shutil.copyfile(found[-1], target)
    return target


def map_file(mtz: Path) -> Path:
    """A 2mFo-DFc map, written as a CCP4 map by gemmi."""
    target = inputs_dir() / "gamma_2fofc.map"
    grid = gemmi.read_mtz_file(str(mtz)).transform_f_phi_to_map("FWT", "PHWT", sample_rate=3)
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(str(target))
    return target


def main():
    mtz = full_mtz()

    i2run("ImportCoordinate", "--XYZIN", f"fullPath={GAMMA / 'gamma_model.pdb'}")
    i2run("ImportSequence", "--SEQIN", f"fullPath={GAMMA / 'gamma.pir'}")
    i2run("ImportAsuContent", "--ASUIN", f"fullPath={GAMMA / 'gamma.asu.xml'}")
    i2run("ImportDictionary", "--DICTIN", f"fullPath={DEMO / 'baz2b' / 'BAZ2BA-x839-LIG.cif'}")
    i2run("ImportUnmerged", "--UNMERGEDIN", f"fullPath={GAMMA / 'HKLOUT_unmerged.mtz'}")

    # One full MTZ, four kinds of data. The columns are named: the file has
    # two sets of map coefficients (FWT/PHWT, DELFWT/PHDELWT), and naming
    # them throughout shows the field in use.
    i2run("ImportObs", "--HKLIN", f"fullPath={mtz}", "--COLUMNS", "F,SIGF")
    i2run("ImportFreeR", "--HKLIN", f"fullPath={mtz}", "--COLUMNS", "FREER")
    i2run("ImportMapCoeffs", "--HKLIN", f"fullPath={mtz}", "--COLUMNS", "FWT,PHWT")
    i2run("ImportPhases", "--HKLIN", f"fullPath={mtz}",
          "--COLUMNS", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB")
    i2run("ImportMap", "--MAPIN", f"fullPath={map_file(mtz)}", "--MAP_SUBTYPE", "1")

    for task in ("ImportCoordinate", "ImportSequence", "ImportAsuContent",
                 "ImportDictionary", "ImportUnmerged", "ImportObs", "ImportFreeR",
                 "ImportMapCoeffs", "ImportPhases", "ImportMap"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
