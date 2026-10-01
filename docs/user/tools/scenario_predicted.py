"""Add to the MDM2 project the run the Process Predicted Models (editbfac)
page is illustrated from.

The crystal holds MDM2's N-terminal domain (97 residues of the 4HG7
construct). The AlphaFold model of human MDM2 (UniProt Q00987) is the whole
491-residue protein, most of it predicted with low confidence: as a search
model it is mostly noise. The question each task answers is the one a user
of a predicted model has: which part of it can be trusted, and does that
part solve the structure?

- editbfac turns pLDDT into B-factors, removes the low-confidence residues
  and splits what is left into compact regions, using the predicted aligned
  error (PAE) where it is given.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_predicted.py

Run scenario_refine.py first (it makes MDM2). The model and its PAE come
from the AlphaFold Database; nothing is sent anywhere. The job gets an
unrun clone.

SliceNDice is not run here. Searching the MDM2 data with this model takes
Phaser over an hour, most of it spent on domains the crystal does not
contain; SliceNDice wants a case of its own (a two-lobed kinase domain).
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

    clone_last(PROJECT, "editbfac")


if __name__ == "__main__":
    main()
