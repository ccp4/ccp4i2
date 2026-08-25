"""
A coordinate file's format must be decided by its content, and each of the two
selection APIs must keep the promise its name makes.

Background in docs/coordinate-format-fidelity.md. The motivating defect:
SubstituteLigand copied an mmCIF input verbatim into a file it named
selected_atoms.pdb, and dimple's reader --- which detects format from the
extension --- failed with "Incorrect file format (perhaps it is cif not pdb?)".
Every assertion here checks what a file *contains*, never what it is called.

Needs gemmi (available CCP4-free), no CCP4 binaries.
"""
import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.core.CCP4ModelData import (
    COORDINATE_FORMAT_MMCIF,
    COORDINATE_FORMAT_PDB,
    CPdbDataFile,
    detect_coordinate_format,
)


def _write(tmp_path, name, fmt, chains=("A", "B")):
    """A small structure, written in *fmt*, under any name we like."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(40, 40, 40, 90, 90, 90)
    model = gemmi.Model("1")
    for chain_name in chains:
        chain = gemmi.Chain(chain_name)
        for seq in (1, 2):
            residue = gemmi.Residue()
            residue.name = "GLY"
            residue.seqid = gemmi.SeqId(seq, " ")
            atom = gemmi.Atom()
            atom.name = "CA"
            atom.element = gemmi.Element("C")
            atom.pos = gemmi.Position(seq, seq, seq)
            residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()

    path = tmp_path / name
    if fmt == COORDINATE_FORMAT_MMCIF:
        st.make_mmcif_document().write_file(str(path))
    else:
        st.write_pdb(str(path))
    return str(path)


def _actual_format(path):
    """What the file really is, read back with an explicit format request."""
    st = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Detect)
    return (COORDINATE_FORMAT_MMCIF if st.input_format == gemmi.CoorFormat.Mmcif
            else COORDINATE_FORMAT_PDB)


def _file(path, selection=None):
    obj = CPdbDataFile()
    obj.setFullPath(path)
    if selection is not None:
        obj.selection.text.set(selection)
    return obj


# --- detection --------------------------------------------------------------

@pytest.mark.parametrize("name,fmt", [
    ("honest.pdb", COORDINATE_FORMAT_PDB),
    ("honest.cif", COORDINATE_FORMAT_MMCIF),
    ("lying.pdb", COORDINATE_FORMAT_MMCIF),   # mmCIF under a .pdb name
    ("lying.cif", COORDINATE_FORMAT_PDB),     # PDB under a .cif name
    ("nameless.tmp", COORDINATE_FORMAT_PDB),  # an extension gemmi does not know
    ("nameless2.tmp", COORDINATE_FORMAT_MMCIF),
])
def test_format_is_detected_from_content(tmp_path, name, fmt):
    path = _write(tmp_path, name, fmt)
    assert detect_coordinate_format(path) == fmt


def test_detection_survives_leading_comments(tmp_path):
    """The hand-rolled sniff read two lines and gave up; mmCIF may open with
    comments or blank lines."""
    path = tmp_path / "commented.cif"
    body = open(_write(tmp_path, "src.cif", COORDINATE_FORMAT_MMCIF)).read()
    path.write_text("# written by something\n#\n\n" + body)
    assert detect_coordinate_format(str(path)) == COORDINATE_FORMAT_MMCIF


def test_an_unreadable_path_is_none(tmp_path):
    assert detect_coordinate_format(str(tmp_path / "absent.pdb")) is None


def test_isMMCIF_agrees_with_the_content_not_the_name(tmp_path):
    assert _file(_write(tmp_path, "lying.pdb", COORDINATE_FORMAT_MMCIF)).isMMCIF() is True
    assert _file(_write(tmp_path, "lying.cif", COORDINATE_FORMAT_PDB)).isMMCIF() is False


# --- getSelectedAtomsPdbFile: PDB, always -----------------------------------

@pytest.mark.parametrize("name,fmt", [
    ("in.pdb", COORDINATE_FORMAT_PDB),
    ("in.cif", COORDINATE_FORMAT_MMCIF),
    ("lying.pdb", COORDINATE_FORMAT_MMCIF),
])
@pytest.mark.parametrize("selection", [None, "A/"])
def test_pdb_api_always_writes_pdb(tmp_path, name, fmt, selection):
    """The case that started this: mmCIF in, no selection, .pdb out."""
    src = _write(tmp_path, name, fmt)
    out = str(tmp_path / "out.pdb")
    assert _file(src, selection).getSelectedAtomsPdbFile(out) == 0
    assert _actual_format(out) == COORDINATE_FORMAT_PDB


def test_pdb_api_keeps_the_selection(tmp_path):
    src = _write(tmp_path, "in.cif", COORDINATE_FORMAT_MMCIF, chains=("A", "B"))
    out = str(tmp_path / "out.pdb")
    assert _file(src, "A/").getSelectedAtomsPdbFile(out) == 0
    st = gemmi.read_structure(out, format=gemmi.CoorFormat.Detect)
    assert [c.name for c in st[0]] == ["A"]


def test_pdb_api_refuses_what_pdb_cannot_express(tmp_path):
    """A chain name PDB has no room for. A refusal beats a truncated model."""
    src = _write(tmp_path, "wide.cif", COORDINATE_FORMAT_MMCIF, chains=("ABCD",))
    out = str(tmp_path / "out.pdb")
    assert _file(src).getSelectedAtomsPdbFile(out) != 0


# --- getSelectedAtomsFile: the input's format, always -----------------------

@pytest.mark.parametrize("name,fmt,expected_ext", [
    ("in.pdb", COORDINATE_FORMAT_PDB, ".pdb"),
    ("in.cif", COORDINATE_FORMAT_MMCIF, ".cif"),
    ("lying.pdb", COORDINATE_FORMAT_MMCIF, ".cif"),
    ("lying.cif", COORDINATE_FORMAT_PDB, ".pdb"),
])
@pytest.mark.parametrize("selection", [None, "A/"])
def test_preserving_api_keeps_the_input_format(tmp_path, name, fmt, expected_ext, selection):
    src = _write(tmp_path, name, fmt)
    out = _file(src, selection).getSelectedAtomsFile("model", str(tmp_path))
    assert out.endswith(expected_ext), "the name must match what was written"
    assert _actual_format(out) == fmt, "and the content must match the input"


def test_the_two_apis_disagree_on_purpose(tmp_path):
    """Same mmCIF input, same selection: one converts, one preserves."""
    src = _write(tmp_path, "in.cif", COORDINATE_FORMAT_MMCIF)
    as_pdb = str(tmp_path / "as.pdb")
    _file(src, "A/").getSelectedAtomsPdbFile(as_pdb)
    preserved = _file(src, "A/").getSelectedAtomsFile("kept", str(tmp_path))
    assert _actual_format(as_pdb) == COORDINATE_FORMAT_PDB
    assert _actual_format(preserved) == COORDINATE_FORMAT_MMCIF


def test_a_missing_input_is_reported_not_copied(tmp_path):
    obj = CPdbDataFile()
    obj.setFullPath(str(tmp_path / "absent.cif"))
    assert obj.getSelectedAtomsPdbFile(str(tmp_path / "out.pdb")) != 0
