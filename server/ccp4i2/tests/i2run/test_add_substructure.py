"""add_substructure: sites put into a model's frame and merged, through i2run."""
import xml.etree.ElementTree as ET

import gemmi

from .utils import demoData, i2run


def _sites_in_another_frame(tmp_path):
    # The sulphurs of 1buh (P21) as an S-SAD substructure, inverted and moved
    # half a cell along c and 37% along the polar b axis
    model = gemmi.read_structure(demoData("CDK1CyclinBCKS2", "1buh.pdb"))
    st = gemmi.Structure()
    st.cell, st.spacegroup_hm = model.cell, model.spacegroup_hm
    chain = gemmi.Chain("S")
    for atom in (a for ch in model[0] for r in ch for a in r if a.element.name == "S"):
        residue = gemmi.Residue()
        residue.name, residue.seqid = "S", gemmi.SeqId(len(chain) + 1, " ")
        site = atom.clone()
        f = model.cell.fractionalize(atom.pos)
        site.pos = model.cell.orthogonalize(gemmi.Fractional(-f.x, -f.y + 0.37, -f.z + 0.5))
        residue.add_atom(site)
        chain.add_residue(residue)
    m = gemmi.Model("1")
    m.add_chain(chain)
    st.add_model(m)
    path = tmp_path / "sites.pdb"
    st.write_pdb(str(path))
    return str(path)


def test_s_sad_sites_found_in_their_frame_and_not_doubled(tmp_path):
    args = ["add_substructure", "--XYZIN", demoData("CDK1CyclinBCKS2", "1buh.pdb"),
            "--XYZIN_SUB", _sites_in_another_frame(tmp_path)]
    with i2run(args) as job:
        done = ET.parse(job / "program.xml").find(".//CompleteModel")
        assert done.find("Origin").get("moved") == "True"
        assert done.get("on_model") == "11" and done.get("added") == "0"
        out = gemmi.read_structure(str(next(job.glob("XYZOUT*.pdb"))))
        model = gemmi.read_structure(demoData("CDK1CyclinBCKS2", "1buh.pdb"))
        assert out[0].count_atom_sites() == model[0].count_atom_sites()
