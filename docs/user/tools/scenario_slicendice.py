"""Build the project the SliceNDice page is illustrated from.

The kinase domain of the Src-family kinase Lck (PDB 4c3f, residues 237-501,
1.72 A, P212121, data and free set from PDB-REDO) solved with the
AlphaFold model of Lck (UniProt P06239). A kinase domain is the case
SliceNDice is for: two lobes joined by a hinge, whose relative orientation a
prediction may not share with the crystal, so the model is tried whole and
split into rigid pieces.

The AlphaFold model is the whole protein (SH3, SH2 and kinase domains); it
is trimmed here to the kinase domain (residues 225-509), as a user would
before searching data that hold only that domain. SliceNDice is asked for
one number of splits: 0.1.3 (CCP4 9) runs MR on only one of the splits it
makes when given a range.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_slicendice.py

Makes the project LckKinase. Nothing is sent anywhere; SliceNDice's own
PDB and AlphaFold DB searches are off. The SliceNDice job gets an unrun
clone.
"""
import json
import urllib.request
import xml.etree.ElementTree as ET

import gemmi

from scenario_common import clone_last, fetch, i2run as _i2run, inputs_dir, scratch_home

PROJECT = "LckKinase"
PDBE = "https://www.ebi.ac.uk/pdbe"
REDO = "https://pdb-redo.eu/db"
KINASE = (225, 509)


def i2run(*args):
    _i2run(PROJECT, *args)


def job_dir(number):
    return scratch_home() / "projects" / PROJECT.lower() / "CCP4_JOBS" / f"job_{number}"


def kinase_domain_model():
    """The AlphaFold model of Lck, kinase domain only."""
    with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/P06239") as r:
        entry = json.load(r)[0]
    full = fetch(entry["pdbUrl"])
    out = inputs_dir() / (full.stem + "_kinase.pdb")
    if not out.exists():
        structure = gemmi.read_structure(str(full))
        for chain in structure[0]:
            for i in reversed(range(len(chain))):
                if not KINASE[0] <= chain[i].seqid.num <= KINASE[1]:
                    del chain[i]
        structure.write_pdb(str(out))
    return out


def main():
    mtz = fetch(f"{REDO}/4c3f/4c3f_final.mtz")
    fasta = fetch(f"{PDBE}/api/v2/pdb/entry/4c3f/fasta", "4c3f.fasta")
    sequence = "".join(fasta.read_text().splitlines()[1:])
    model = kinase_domain_model()

    i2run("import_merged", "--HKLIN", f"fullPath={mtz}",
          "--HKLIN_OBS_COLUMNS", "FP,SIGFP", "--HKLIN_FREER_COLUMN", "FREE")
    i2run("ProvideAsuContents",
          "--ASU_CONTENT", f"sequence={sequence}", "nCopies=1",
          "name=LCK", "description=Lck kinase domain (PDB 4c3f)",
          "polymerType=PROTEIN", "source/baseName=4c3f.fasta",
          f"source/relPath={fasta.parent}")

    # Two splits: the N- and C-lobes, each a search model.
    i2run("slicendice",
          "--F_SIGF", "fileOut=import_merged[-1].OBSOUT",
          "--FREERFLAG", "fileOut=import_merged[-1].FREEOUT",
          "--ASUIN", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE",
          "--XYZIN", f"fullPath={model}",
          "--BFACTOR_TREATMENT", "plddt",
          "--SEARCH_PDB", "False", "--SEARCH_AFDB", "False",
          "--NO_MOLS", "1", "--MIN_SPLITS", "2", "--MAX_SPLITS", "2")
    best = ET.parse(job_dir(3) / "program.xml").find(".//RunInfo/Best")
    assert best is not None and best.findtext("Solved") == "True", \
        f"SliceNDice did not solve it: {ET.tostring(best) if best is not None else 'no results'}"

    clone_last(PROJECT, "slicendice")


if __name__ == "__main__":
    main()
