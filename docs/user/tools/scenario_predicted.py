"""Add to the MDM2 project the runs the predicted-model pages are illustrated
from: Process Predicted Models (editbfac) and SliceNDice.

The crystal holds MDM2's N-terminal domain (97 residues of the 4HG7
construct). The AlphaFold model of human MDM2 (UniProt Q00987) is the whole
491-residue protein, most of it predicted with low confidence: as a search
model it is mostly noise. The question each task answers is the one a user
of a predicted model has: which part of it can be trusted, and does that
part solve the structure?

- editbfac turns pLDDT into B-factors, removes the low-confidence residues
  and splits what is left into compact regions, using the predicted aligned
  error (PAE) where it is given.
- SliceNDice does the same preparation itself and runs Phaser with each
  split, here against the MDM2 data (1.35 A, P6522).

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_predicted.py

Run scenario_refine.py first (it makes MDM2: data in job 1, AU contents in
job 9). The model and its PAE come from the AlphaFold Database; nothing is
sent anywhere. SliceNDice's own searches of the PDB and the AlphaFold DB
are turned off. Each task's last job gets an unrun clone.
"""
import json
import urllib.request

from scenario_common import clone_last, fetch, i2run

PROJECT = "MDM2"
UNIPROT = "Q00987"


def alphafold():
    """The model and PAE file of the current AlphaFold DB entry."""
    with urllib.request.urlopen(
            f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}") as response:
        entry = json.load(response)[0]
    return fetch(entry["pdbUrl"]), fetch(entry["paeDocUrl"])


def main():
    model, pae = alphafold()

    i2run(PROJECT, "editbfac", "--XYZIN", f"fullPath={model}",
          "--PAEIN", f"fullPath={pae}")

    # Job numbers in MDM2: 1 the data reduction (data and free set), 9 the
    # AU contents. One copy of the domain in the asymmetric unit.
    i2run(PROJECT, "slicendice",
          "--F_SIGF", "fileOut=[1].HKLOUT[0]",
          "--FREERFLAG", "fileOut=[1].FREEROUT",
          "--ASUIN", "fileOut=[9].ASUCONTENTFILE",
          "--XYZIN", f"fullPath={model}",
          "--BFACTOR_TREATMENT", "plddt",
          "--SEARCH_PDB", "False", "--SEARCH_AFDB", "False",
          "--NO_MOLS", "1")

    for task in ("editbfac", "slicendice"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
