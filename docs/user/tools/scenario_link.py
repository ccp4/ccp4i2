"""Build the project the Make Covalent Link page is illustrated from.

PDB entry 6ndn: pyridoxal 5'-phosphate (PLP) bound to a lysine as a Schiff
base, the internal aldimine of a PLP enzyme. The link joins the lysine's NZ
to PLP's C4A by a double bond; PLP loses the aldehyde oxygen O4A in the
condensation.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_link.py

Makes the project PLP_Link: one MakeLink run, and an unrun clone.
"""
import gemmi

from scenario_common import clone_last, fetch, i2run, scratch_home

PROJECT = "PLP_Link"


def main():
    model = fetch("https://www.ebi.ac.uk/pdbe/entry-files/download/pdb6ndn.ent",
                  "6ndn.pdb")
    i2run(PROJECT, "MakeLink",
          "--RES_NAME_1_TLC", "LYS", "--RES_NAME_2_TLC", "PLP",
          "--ATOM_NAME_1_TLC", "NZ", "--ATOM_NAME_2_TLC", "C4A",
          "--ATOM_NAME_1", "NZ", "--ATOM_NAME_2", "C4A",
          "--TOGGLE_DELETE_2", "True", "--DELETE_2", "O4A",
          "--BOND_ORDER", "DOUBLE",
          "--TOGGLE_LINK", "True", "--XYZIN", f"fullPath={model}")
    # The story: the model comes back carrying the link.
    job = scratch_home() / "projects/plp_link/CCP4_JOBS/job_1"
    structure = gemmi.read_pdb(str(job / "ModelWithLinks.pdb"))
    assert any(c.link_id == "LYS-PLP" for c in structure.connections), \
        "no LYS-PLP link in the output model"
    clone_last(PROJECT, "MakeLink")


if __name__ == "__main__":
    main()
