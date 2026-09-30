"""Build the project the SubstituteLigand help page is illustrated from.

MDM2 with Nutlin-3a, from the demo data that ships with CCP4i2: the parent
structure 4hg7 with its ligand and waters removed by the atom selection, the
unmerged reflections of a soak, and Nutlin-3a given as a SMILES string.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_SubstituteLigand.py

Makes the project MDM2_Nutlin: one SubstituteLigand run, and an unrun clone.
"""
from pathlib import Path

import gemmi

from scenario_common import clone_last, i2run, scratch_home

PROJECT = "MDM2_Nutlin"
DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data/mdm2"
# Nutlin-3a, with its stereochemistry (demo_data/mdm2/README.html).
NUTLIN_3A = ("COc1ccc(c(OC(C)C)c1)C2=N[C@H]([C@H](N2C(=O)N3CCNC(=O)C3)"
             "c4ccc(Cl)cc4)c5ccc(Cl)cc5")


def main():
    i2run(PROJECT, "SubstituteLigand",
          "--XYZIN", f"fullPath={DEMO / '4hg7.cif'}",
          "selection/text=not (NUT) and not (HOH)",
          "--UNMERGEDFILES", f"file={DEMO / 'mdm2_unmerged.mtz'}",
          "--SMILESIN", NUTLIN_3A,
          "--PIPELINE", "DIMPLE")
    # The page's story is that the ligand is found: check it, rather than
    # illustrate a run that quietly placed nothing (an atom selection that
    # matched nothing once left Nutlin in the site, and DRG went unplaced).
    out = scratch_home() / "projects/mdm2_nutlin/CCP4_JOBS/job_1/XYZOUT.pdb"
    names = {r.name for ch in gemmi.read_structure(str(out))[0] for r in ch}
    assert "DRG" in names and "NUT" not in names, sorted(names)
    clone_last(PROJECT, "SubstituteLigand")


if __name__ == "__main__":
    main()
