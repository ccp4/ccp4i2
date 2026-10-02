"""Build the runs the remaining refinement and validation pages are illustrated
from.

- lorestr_i2 (CDK2_CyclinA): refinement of the CDK2-only partial model of
  1h1s with external restraints from a local reference, 1jst (CDK2/cyclin A
  from another crystal, in demo_data), fetching nothing (AUTO none).
- dnatco_pipe and nucleofind (RNA_1hr2, a new project): the P4-P6 RNA domain
  1hr2. DNATCO classifies the conformation of each dinucleotide in the
  deposited model and compares it with the PDB-REDO re-refinement; NucleoFind
  predicts phosphate, sugar and base positions from the PDB-REDO map.
- metalCoord (MetalSite, a new project): the AlF3 of 4dl8, restraints from
  the metal-coordination statistics of the PDB, at the task's defaults.
- buster and pdb_redo_api (CDK2_CyclinA): set up and not run (i2run --delay).
  BUSTER is not installed here; PDB-REDO is a web service, and nothing is
  sent to it.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_refine_tail.py

Needs scenario_parrot.py (CDK2_CyclinA) first. Downloads 1hr2 (RCSB,
PDB-REDO) and 4dl8 (PDBe) once into the scratch home. Each run gets an unrun
clone. Name tasks on the command line to run only those.
"""
import sys
from pathlib import Path

import ccp4i2
from scenario_common import clone_last, fetch, i2run, last_job_dir

DEMO = Path(ccp4i2.__file__).parent / "demo_data" / "CDK1CyclinBCKS2"
RCSB = "https://files.rcsb.org/download"
REDO = "https://pdb-redo.eu/db"
PDBE = "https://www.ebi.ac.uk/pdbe/entry-files/download"

CDK2 = ["--F_SIGF", "fileIn=[2].F_SIGF", "--FREERFLAG", "fileIn=[2].FREERFLAG",
        "--XYZIN", "fileOut=[2].XYZOUT"]


def runs():
    rna = fetch(f"{RCSB}/1hr2.cif", "1hr2.cif")
    rna_redo = fetch(f"{REDO}/1hr2/1hr2_final.pdb", "1hr2_final.pdb")
    rna_mtz = fetch(f"{REDO}/1hr2/1hr2_final.mtz", "1hr2_final.mtz")
    metal = fetch(f"{PDBE}/4dl8.cif", "4dl8.cif")
    return [
        ("CDK2_CyclinA", "lorestr_i2", CDK2 + [
            "--REFERENCE_LIST", f"fullPath={DEMO / '1jst.pdb'}",
            "--AUTO", "none", "--CPU", "4"], True),
        ("RNA_1hr2", "dnatco_pipe", [
            "--XYZIN1", f"fullPath={rna}", "--TOGGLE_XYZIN2", "True",
            "--XYZIN2", f"fullPath={rna_redo}"], True),
        ("RNA_1hr2", "nucleofind", [
            "--FPHIIN", f"fullPath={rna_mtz}", "columnLabels=/*/*/[FWT,PHWT]"], True),
        ("MetalSite", "metalCoord", [
            "--XYZIN", f"fullPath={metal}", "--LIGAND_CODE", "AF3"], True),
        ("CDK2_CyclinA", "buster", ["--delay"] + CDK2, False),
        ("CDK2_CyclinA", "pdb_redo_api", ["--delay"] + CDK2, False),
    ]


def main(only=()):
    failed = []
    for project, task, args, run in runs():
        if only and task not in only:
            continue
        try:
            i2run(project, task, *args)
        except Exception as err:  # keep going: one task's failure is a finding
            failed.append(f"{project} {task}: {err}")
            continue
        print(task, "->", last_job_dir(project, task))
        if run:
            clone_last(project, task)
    if failed:
        raise SystemExit("Failed:\n" + "\n".join(failed))


if __name__ == "__main__":
    main(sys.argv[1:])
