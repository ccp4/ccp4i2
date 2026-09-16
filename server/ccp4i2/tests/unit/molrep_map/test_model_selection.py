"""molrep_map honours an atom selection on the model it places.

So a user can fetch a reference structure (e.g. 3g33) and place only the chains
they want (e.g. A and D) via the task interface. _ensure_pdb routes the model
through getSelectedAtomsPdbFile, which applies any selection *and* gives molrep
the PDB it needs -- whole model, converted, when nothing is selected.

Needs gemmi (CCP4-free), no molrep binary.
"""

import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.core.tasks import get_plugin_class


def _model(tmp_path, name, chains, fmt="pdb"):
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(40, 40, 40, 90, 90, 90)
    model = gemmi.Model("1")
    for cname in chains:
        chain = gemmi.Chain(cname)
        for seq in (1, 2, 3):
            res = gemmi.Residue()
            res.name = "GLY"
            res.seqid = gemmi.SeqId(seq, " ")
            atom = gemmi.Atom()
            atom.name = "CA"
            atom.element = gemmi.Element("C")
            atom.pos = gemmi.Position(seq, seq, seq)
            res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    path = tmp_path / name
    if fmt == "cif":
        st.make_mmcif_document().write_file(str(path))
    else:
        st.write_pdb(str(path))
    return str(path)


def _ensure(tmp_path, src, selection=None):
    p = get_plugin_class("molrep_map")()
    p.workDir = str(tmp_path)
    xyzin = p.container.inputData.XYZIN
    xyzin.setFullPath(src)
    if selection is not None:
        xyzin.selection.text.set(selection)
    out = p._ensure_pdb(xyzin)
    return out, gemmi.read_structure(out)


def _chains(structure):
    return sorted(ch.name for ch in structure[0])


def test_selection_keeps_only_the_chosen_chains(tmp_path):
    # The motivating case: fetch a reference (chains A, B, D) and place only A+D.
    src = _model(tmp_path, "ref.pdb", ["A", "B", "D"])
    out, st = _ensure(tmp_path, src, selection="A,D/")
    assert _chains(st) == ["A", "D"]
    assert out.endswith(".pdb")


def test_single_chain_selection(tmp_path):
    src = _model(tmp_path, "ref.pdb", ["A", "B", "D"])
    _, st = _ensure(tmp_path, src, selection="A/")
    assert _chains(st) == ["A"]


def test_no_selection_keeps_every_chain(tmp_path):
    src = _model(tmp_path, "ref.pdb", ["A", "B", "D"])
    _, st = _ensure(tmp_path, src)
    assert _chains(st) == ["A", "B", "D"]


def test_mmcif_input_is_converted_to_pdb(tmp_path):
    # No selection, mmCIF in -> molrep still gets PDB out.
    src = _model(tmp_path, "ref.cif", ["A", "B"], fmt="cif")
    out, st = _ensure(tmp_path, src)
    assert out.endswith(".pdb")
    assert gemmi.read_structure(out, format=gemmi.CoorFormat.Detect).input_format \
        == gemmi.CoorFormat.Pdb
    assert _chains(st) == ["A", "B"]
