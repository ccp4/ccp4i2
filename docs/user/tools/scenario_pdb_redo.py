"""Build the run the PDB-REDO page is illustrated from.

- pdb_redo_api (PDB_REDO_2d2k, a new project): the minimal hairpin ribozyme
  2d2k (2.65 A, P6122), re-refined and rebuilt by the PDB-REDO web service.
  One of the cases Robbie Joosten (PDB-REDO) suggested, and the one among
  them that deposited its free set, so the run can reuse the free set the
  model was refined against, as the page tells the reader to.

Unlike every other scenario, this one sends data to an outside service:
a deposited entry's model, data and sequence, to pdb-redo.eu, with the
PDB-REDO token of whoever runs it (Preferences > Credentials, or
PDB_REDO_TOKEN_ID / PDB_REDO_TOKEN_SECRET). Only ever run it on public data.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_pdb_redo.py

The model goes in as mmCIF, never PDB format: PDB-REDO is dropping
PDB-format output. Downloads 2d2k (RCSB, PDBe) once into the scratch home.
A run takes tens of minutes at the service. Ends with an unrun clone.
"""
import json

from scenario_common import clone_last, fetch, i2run, last_job_dir

PROJECT = "PDB_REDO_2d2k"
RCSB = "https://files.rcsb.org/download"
PDBE = "https://www.ebi.ac.uk/pdbe"


def main():
    model = fetch(f"{RCSB}/2d2k.cif", "2d2k.cif")
    sf = fetch(f"{RCSB}/2d2k-sf.cif", "2d2k-sf.cif")
    fasta = fetch(f"{PDBE}/api/v2/pdb/entry/2d2k/fasta", "2d2k.fasta")

    # 1. The deposited data, with the deposited free set (status 'f').
    i2run(PROJECT, "import_merged", "--HKLIN", f"fullPath={sf}")

    # 2. PDB-REDO on the deposited model, at the task's defaults.
    i2run(PROJECT, "pdb_redo_api",
          "--XYZIN", f"fullPath={model}",
          "--F_SIGF", "fileOut=import_merged[-1].OBSOUT",
          "--FREERFLAG", "fileOut=import_merged[-1].FREEOUT",
          "--SEQIN", f"fullPath={fasta}")

    # The run must have brought back both models, as mmCIF, and its R factors.
    # PDB-REDO names its files after its own run code (Kjeb_final.cif), not
    # the entry, and this run wrote PDB-format copies too: the mmCIF is taken.
    job = last_job_dir(PROJECT, "pdb_redo_api")
    for model in ("final", "besttls"):
        assert list(job.glob(f"*_{model}.cif")), f"no {model} mmCIF in {job}"
        assert not list(job.glob(f"*_{model}.pdb")), f"{model} taken as PDB format"
    data = json.loads(next(job.rglob("data.json")).read_text())["properties"]
    print("PDB-REDO R/R-free: deposited", data["RFACT"], data["RFREE"],
          "recalculated", data["RCAL"], data["RFCAL"],
          "final", data["RFIN"], data["RFFIN"])
    assert data["RFFIN"] < data["RFCAL"], "R-free did not improve"
    clone_last(PROJECT, "pdb_redo_api")


if __name__ == "__main__":
    main()
