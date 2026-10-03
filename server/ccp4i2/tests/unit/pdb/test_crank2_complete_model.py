"""Crank2's complete model: the built model plus the substructure, once each.

Refining Crank2's XYZOUT without its substructure lost 0.15 in R-free for a
mercury soak (HypF). A site that is already a model atom must not be added
again: S-SAD sites are cysteine and methionine sulphurs, and a SeMet site is
the Se of a residue built as methionine.
"""
import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.pipelines.crank2.script.complete_model import complete_model

CELL = gemmi.UnitCell(40, 50, 60, 90, 90, 90)


def _structure(atoms, chain="A"):
    """atoms: (residue name, seqnum, atom name, element, (x, y, z), occupancy)"""
    st = gemmi.Structure()
    st.cell = CELL
    st.spacegroup_hm = "P 21 21 21"
    model = gemmi.Model("1")
    ch = gemmi.Chain(chain)
    residues = {}
    for resname, seqnum, name, element, xyz, occ in atoms:
        key = (resname, seqnum)
        if key not in residues:
            residue = gemmi.Residue()
            residue.name = resname
            residue.seqid = gemmi.SeqId(seqnum, " ")
            residues[key] = residue
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element(element)
        atom.pos = gemmi.Position(*xyz)
        atom.occ = occ
        atom.b_iso = 20.0
        residues[key].add_atom(atom)
    for residue in residues.values():
        ch.add_residue(residue)
    model.add_chain(ch)
    st.add_model(model)
    st.setup_entities()
    return st


MODEL = [
    ("CYS", 1, "CA", "C", (10.0, 10.0, 10.0), 1.0),
    ("CYS", 1, "SG", "S", (11.5, 10.0, 10.0), 1.0),
    ("MET", 2, "CA", "C", (20.0, 20.0, 20.0), 1.0),
    ("MET", 2, "SD", "S", (21.8, 20.0, 20.0), 1.0),
    ("ALA", 3, "CA", "C", (30.0, 30.0, 30.0), 1.0),
]


def _run(tmp_path, sites):
    model = tmp_path / "model.pdb"
    substr = tmp_path / "substr.pdb"
    out = tmp_path / "complete.pdb"
    _structure(MODEL).write_pdb(str(model))
    _structure(sites, chain="S").write_pdb(str(substr))
    report = complete_model(model, substr, out)
    return report, gemmi.read_structure(str(out))[0]


def _atoms(model):
    return [(ch.name, r.name, a.name, a.element.name) for ch in model for r in ch for a in r]


def test_soak_sites_are_added_with_their_occupancy(tmp_path):
    report, model = _run(tmp_path, [("HG", 1, "HG", "Hg", (5.0, 30.0, 45.0), 0.69),
                                   ("HG", 2, "HG", "Hg", (35.0, 5.0, 5.0), 0.46)])
    assert report["added"] == ["HG occupancy 0.69", "HG occupancy 0.46"]
    hg = [a for ch in model for r in ch for a in r if a.element.name == "Hg"]
    assert [round(a.occ, 2) for a in hg] == [0.69, 0.46]
    assert len(_atoms(model)) == len(MODEL) + 2


def test_sulphur_sites_on_the_model_are_not_added_again(tmp_path):
    report, model = _run(tmp_path, [("S", 1, "S", "S", (11.6, 10.1, 10.0), 0.9),
                                   ("S", 2, "S", "S", (21.7, 20.0, 20.1), 0.8)])
    assert report["added"] == []
    assert len(report["on_model"]) == 2
    assert len(_atoms(model)) == len(MODEL)


def test_selenium_on_a_methionine_makes_it_selenomethionine(tmp_path):
    report, model = _run(tmp_path, [("SE", 1, "SE", "Se", (21.9, 20.1, 20.0), 0.9)])
    assert report["added"] == []
    assert report["converted"] == ["A/MET2"]
    names = _atoms(model)
    assert ("A", "MSE", "SE", "Se") in names
    assert not any(r == "MET" for _, r, _, _ in names)
    assert len(names) == len(MODEL)


def test_a_site_on_top_of_another_atom_is_left_out_and_said(tmp_path):
    report, model = _run(tmp_path, [("XE", 1, "XE", "Xe", (21.0, 20.0, 20.0), 0.3)])
    assert report["added"] == []
    assert len(report["clashes"]) == 1 and "MET2" in report["clashes"][0]
    assert len(_atoms(model)) == len(MODEL)


def test_a_site_on_a_symmetry_copy_of_a_model_atom_is_recognised(tmp_path):
    # P 21 21 21 operator (-x+1/2, -y, z+1/2) takes SG (11.5, 10, 10) to
    # (8.5, -10, 40); a site there is that cysteine's sulphur
    report, model = _run(tmp_path, [("S", 1, "S", "S", (8.5, -10.0, 40.0), 0.9)])
    assert report["added"] == [] and len(report["on_model"]) == 1


def test_added_sites_go_in_a_chain_of_their_own(tmp_path):
    report, model = _run(tmp_path, [("IOD", 1, "I", "I", (5.0, 30.0, 45.0), 0.5)])
    chains = {ch.name: [r.name for r in ch] for ch in model}
    assert chains["A"] == ["CYS", "MET", "ALA"]
    new = [name for name in chains if name != "A"]
    assert len(new) == 1 and chains[new[0]] == ["IOD"]
