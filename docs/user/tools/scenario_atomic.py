"""Build the project the atomic-resolution phasing pages are illustrated from:
Single Atom Molecular Replacement (phaser_singleMR) and ACORN.

PDB entry 3njw, the first high-resolution structure of a lasso peptide: 19
residues, two of them cysteines, measured to 0.86 A (P212121), with the
data and free set from PDB-REDO. At that resolution the structure can be
solved from almost nothing:

- phaser_singleMR places single atoms as if they were a search model (here
  the two sulfurs) and completes the substructure from log-likelihood-gain
  maps until the model is the whole peptide;
- ACORN starts from those phases and refines and extends them by dynamic
  density modification, the ab initio route for data this good.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_atomic.py

Makes the project LassoPeptide. Each task's last job gets an unrun clone.
"""
import re
import xml.etree.ElementTree as ET

from scenario_common import clone_last, fetch, i2run as _i2run, scratch_home

PROJECT = "LassoPeptide"
PDBE = "https://www.ebi.ac.uk/pdbe"
REDO = "https://pdb-redo.eu/db"


def i2run(*args):
    _i2run(PROJECT, *args)


def job_dir(number):
    return scratch_home() / "projects" / PROJECT.lower() / "CCP4_JOBS" / f"job_{number}"


def main():
    mtz = fetch(f"{REDO}/3njw/3njw_final.mtz")
    fasta = fetch(f"{PDBE}/api/v2/pdb/entry/3njw/fasta", "3njw.fasta")
    sequence = "".join(fasta.read_text().splitlines()[1:])
    assert len(sequence) == 19, sequence

    # 1. The data and PDB-REDO's free set.
    i2run("import_merged", "--HKLIN", f"fullPath={mtz}",
          "--HKLIN_OBS_COLUMNS", "FP,SIGFP", "--HKLIN_FREER_COLUMN", "FREE")

    # 2. AU contents: one peptide.
    i2run("ProvideAsuContents",
          "--ASU_CONTENT", f"sequence={sequence}", "nCopies=1",
          "name=LASSO", "description=Lasso peptide (PDB 3njw)",
          "polymerType=PROTEIN", "source/baseName=3njw.fasta",
          f"source/relPath={fasta.parent}")

    # 3. Two sulfur atoms, then LLG completion to the whole peptide.
    i2run("phaser_singleMR",
          "--F_SIGF", "fileOut=import_merged[-1].OBSOUT",
          "--FREERFLAG", "fileOut=import_merged[-1].FREEOUT",
          "--COMP_BY", "ASU", "--ASUFILE", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE",
          "--SINGLE_ATOM_TYPE", "S", "--SINGLE_ATOM_NUM", "2")
    log = (job_dir(3) / "log.txt").read_text()
    rwork = float(re.findall(r"Final R-factor = +([\d.]+)", log)[-1])
    assert rwork < 25, f"single-atom MR did not solve it: R {rwork}"

    # 4. ACORN from the single-atom solution's phases.
    i2run("acorn",
          "--F_SIGF", "fileOut=import_merged[-1].OBSOUT",
          "--ACORN_PHSIN_TYPE", "phases",
          "--ABCD", "fileOut=phaser_singleMR[-1].ABCDOUT[0]")
    cc = float(ET.parse(job_dir(4) / "program.xml").findall(".//CorrelationCoef")[-1].text)
    assert cc > 0.1, f"ACORN final CC {cc}"

    for task in ("phaser_singleMR", "acorn"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
