"""Add to the MDM2 project the runs the scripted Coot pages are illustrated
from: completing a molecular-replacement model.

Job 25 is Phaser's placement of MDMX pruned by Chainsaw into the MDM2 data:
placed clearly (TFZ 12.7), but a homologue with its side chains truncated,
refined only to R-free 0.50. Each task does one step of completing it:

- coot_rsr_morph moves the placed model locally into its map (real-space
  refinement morphing: each region shifts with its neighbours);
- coot_script_lines, from the "Fill partial residues" starting point, gives
  back the side chains Chainsaw truncated;
- prosmart_refmac refines the result, against the same data and free set;
- coot_find_waters adds waters to the refined model.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_coot.py

Run scenario_refine.py and scenario_mr_steps.py first. Each step asserts
its outcome; the Coot tasks run Coot without graphics. Each Coot task's job
gets an unrun clone.
"""
import re

import gemmi

from scenario_common import clone_last, i2run, last_job_dir

PROJECT = "MDM2"

# The "Fill partial residues" starting point (coot_script_templates.ts).
FILL_PARTIAL = """# Input molecules are given internal identifiers "MolHandle_1" etc
# Input maps are given internal identifiers "MapHandle_1" etc
# Input difference maps are given internal identifiers "DifmapHandle_1" etc
#
# The appropriate place to put output pdb files to be picked up by the system
# is stored in the variable dropDir (see below for how to use this)
#
# Beware spaces...this is python after all
fill_partial_residues(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
"""


def atoms(path):
    structure = gemmi.read_structure(str(path))
    return sum(1 for chain in structure[0] for residue in chain for _ in residue)


def output_model(job_dir):
    files = sorted(job_dir.glob("*.pdb")) + sorted(job_dir.glob("*.cif"))
    assert files, f"no model in {job_dir}"
    return max(files, key=lambda f: f.stat().st_mtime)


def r_free(job_dir):
    values = re.findall(r"<r_free>\s*([\d.]+)", (job_dir / "program.xml").read_text())
    return float(values[-1])


def main():
    placed = "fileOut=[25].XYZOUT[0]"
    phaser_map = "fileOut=[25].MAPOUT[0]"
    before = r_free(last_job_dir(PROJECT, "phaser_simple_phil") / "job_3")

    i2run(PROJECT, "coot_rsr_morph", "--XYZIN", placed, "--FPHIIN", phaser_map)
    morphed = output_model(last_job_dir(PROJECT, "coot_rsr_morph"))

    i2run(PROJECT, "coot_script_lines",
          "--XYZIN", "fileOut=coot_rsr_morph[-1].XYZOUT",
          "--FPHIIN", phaser_map,
          "--STARTPOINT", "FILL_PARTIAL_RESIDUES", "--SCRIPT", FILL_PARTIAL)
    filled = output_model(last_job_dir(PROJECT, "coot_script_lines"))
    assert atoms(filled) > atoms(morphed), \
        f"fill_partial_residues added nothing: {atoms(morphed)} -> {atoms(filled)}"

    i2run(PROJECT, "prosmart_refmac",
          "--F_SIGF", "fileOut=[1].HKLOUT[0]", "--FREERFLAG", "fileOut=[1].FREEROUT",
          "--XYZIN", "fileOut=coot_script_lines[-1].XYZOUT[0]")
    after = r_free(last_job_dir(PROJECT, "prosmart_refmac"))
    assert after < before, f"R-free {before} -> {after}: completing did not help"

    i2run(PROJECT, "coot_find_waters",
          "--XYZIN", "fileOut=prosmart_refmac[-1].XYZOUT",
          "--FPHIIN", "fileOut=prosmart_refmac[-1].FPHIOUT")
    watered = output_model(last_job_dir(PROJECT, "coot_find_waters"))
    waters = sum(1 for c in gemmi.read_structure(str(watered))[0] for r in c if r.name == "HOH")
    assert waters > 0, "find waters found none"
    print(f"R-free {before:.3f} -> {after:.3f}; atoms {atoms(morphed)} -> {atoms(filled)}; "
          f"{waters} waters")

    for task in ("coot_rsr_morph", "coot_script_lines", "coot_find_waters"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
