"""The atom-selection preview and checks on CPdbDataFile (#876).

selectionSummary feeds the preview under the selection builder: atoms and
residues per chain, split into polymer, water and other, so a user who ticks a
chain sees that its waters and ligands come with it. validity() refuses a
selection the parser cannot read, and the run-time check refuses one that
matches no atoms.

Needs gemmi (CCP4-free), no binaries.
"""

import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.core.base_object.error_reporting import SEVERITY_ERROR
from ccp4i2.core.tasks import get_plugin_class


def _residue(name, seq, atoms, het=False):
    res = gemmi.Residue()
    res.name = name
    res.seqid = gemmi.SeqId(seq, " ")
    res.het_flag = "H" if het else "A"
    for i, (aname, el) in enumerate(atoms):
        atom = gemmi.Atom()
        atom.name = aname
        atom.element = gemmi.Element(el)
        atom.pos = gemmi.Position(seq, i, 0)
        res.add_atom(atom)
    return res


def _model(tmp_path):
    """Chain A: two GLY (4 atoms each), an ATP-like ligand (3 atoms), two waters.
    Chain B: three GLY."""
    gly = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(40, 40, 40, 90, 90, 90)
    model = gemmi.Model("1")
    a = gemmi.Chain("A")
    a.add_residue(_residue("GLY", 1, gly))
    a.add_residue(_residue("GLY", 2, gly))
    a.add_residue(_residue("LIG", 101, [("C1", "C"), ("C2", "C"), ("O1", "O")], het=True))
    a.add_residue(_residue("HOH", 201, [("O", "O")], het=True))
    a.add_residue(_residue("HOH", 202, [("O", "O")], het=True))
    model.add_chain(a)
    b = gemmi.Chain("B")
    for seq in (1, 2, 3):
        b.add_residue(_residue("GLY", seq, gly))
    model.add_chain(b)
    st.add_model(model)
    path = tmp_path / "model.pdb"
    st.write_pdb(str(path))
    return str(path)


@pytest.fixture
def plugin(tmp_path):
    p = get_plugin_class("molrep_mr")()
    p.workDir = str(tmp_path)
    p.container.inputData.XYZIN.setFullPath(_model(tmp_path))
    return p


def _errors(report, code):
    return [e for e in report.getErrors() if e.get("code") == code]


def test_summary_splits_a_chain_into_polymer_water_and_other(plugin):
    s = plugin.container.inputData.XYZIN.selectionSummary("A/")
    assert s["status"] == "ok"
    assert s["nAtoms"] == 8 + 3 + 2
    assert (s["nResidues"], s["polymer"], s["water"], s["other"]) == (5, 2, 2, 1)
    assert [c["id"] for c in s["chains"]] == ["A"]


def test_summary_of_everything_lists_each_chain(plugin):
    s = plugin.container.inputData.XYZIN.selectionSummary("")
    assert s["status"] == "ok"
    assert s["nAtoms"] == 13 + 12
    assert {c["id"]: c["nResidues"] for c in s["chains"]} == {"A": 5, "B": 3}


def test_summary_uses_the_files_own_selection_when_given_none(plugin):
    xyzin = plugin.container.inputData.XYZIN
    xyzin.selection.text.set("B/")
    s = xyzin.selectionSummary()
    assert s["status"] == "ok"
    assert s["nAtoms"] == 12


def test_summary_of_a_selection_matching_nothing(plugin):
    s = plugin.container.inputData.XYZIN.selectionSummary("Z/")
    assert s["status"] == "ok"
    assert s["nAtoms"] == 0
    assert s["chains"] == []


def test_summary_reports_a_syntax_error_without_reading_the_file(plugin):
    s = plugin.container.inputData.XYZIN.selectionSummary("{A/ or")
    assert s["status"] == "syntax"
    assert s["error"]


def test_a_list_of_residue_names_parses(plugin):
    """The builder writes "(ALA,GLY)" for its residue-name filter; the parser
    named a token type that does not exist on reaching the comma."""
    from ccp4i2.core.CCP4ModelData import CPdbDataFile
    assert CPdbDataFile.selectionSyntaxError("(GLY,LIG)") is None
    s = plugin.container.inputData.XYZIN.selectionSummary("A and (GLY,LIG)")
    assert s["status"] == "ok"
    assert (s["polymer"], s["other"], s["water"]) == (2, 1, 0)


def test_summary_without_a_file_is_unavailable(tmp_path):
    p = get_plugin_class("molrep_mr")()
    assert p.container.inputData.XYZIN.selectionSummary("A/")["status"] == "unavailable"


def test_validity_refuses_an_unreadable_selection(plugin):
    xyzin = plugin.container.inputData.XYZIN
    xyzin.selection.text.set("{A/ or")
    found = _errors(xyzin.validity(), 420)
    assert len(found) == 1
    assert found[0]["severity"] == SEVERITY_ERROR


def test_validity_accepts_a_readable_selection(plugin):
    xyzin = plugin.container.inputData.XYZIN
    xyzin.selection.text.set("A/")
    assert _errors(xyzin.validity(), 420) == []


def test_validity_ignores_a_selection_the_task_does_not_offer(tmp_path):
    """Where ifAtomSelection is off the user cannot see the selection to fix it."""
    from ccp4i2.core.CCP4ModelData import CPdbDataFile
    f = CPdbDataFile(name="XYZIN")
    f.selection.text.set("{A/ or")
    assert _errors(f.validity(), 420) == []


def test_run_time_validity_refuses_a_selection_matching_no_atoms(plugin):
    plugin.container.inputData.XYZIN.selection.text.set("Z/")
    found = _errors(plugin.runTimeValidity(), 221)
    assert len(found) == 1
    assert found[0]["severity"] == SEVERITY_ERROR


def test_run_time_validity_accepts_a_selection_matching_atoms(plugin):
    plugin.container.inputData.XYZIN.selection.text.set("A/")
    assert _errors(plugin.runTimeValidity(), 221) == []
