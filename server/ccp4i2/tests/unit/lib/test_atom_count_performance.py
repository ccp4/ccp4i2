"""
CAtomCountPerformance.setFromPdbDataFile -- the KPI chainsaw and sculptor
have always asked for, and which was never ported from the Qt code.

Both wrappers called it from processOutputFiles(); both raised
``AttributeError: 'CAtomCountPerformance' object has no attribute
'setFromPdbDataFile'`` on every run, and both jobs were recorded as Finished.
C1 of docs/error-handling-remediation.md made the exception visible.

Needs gemmi (available CCP4-free), no CCP4 binaries.
"""
import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.core.CCP4PerformanceData import CAtomCountPerformance


def _structure(tmp_path, name, fmt):
    """A two-residue, three-atom structure written as PDB or mmCIF."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for seq, (resname, atoms) in enumerate(
        [("GLY", ["N", "CA"]), ("HOH", ["O"])], start=1
    ):
        residue = gemmi.Residue()
        residue.name = resname
        residue.seqid = gemmi.SeqId(seq, " ")
        for atom_name in atoms:
            atom = gemmi.Atom()
            atom.name = atom_name
            atom.element = gemmi.Element("C" if atom_name != "O" else "O")
            atom.pos = gemmi.Position(seq, seq, seq)
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()

    path = tmp_path / name
    if fmt == "mmcif":
        st.make_mmcif_document().write_file(str(path))
    else:
        st.write_pdb(str(path))
    return str(path)


class _FileStub:
    """Stands in for a CPdbDataFile: what the wrappers actually pass."""
    def __init__(self, path):
        self._path = path

    def getFullPath(self):
        return self._path

    def __str__(self):
        return self._path


def test_counts_atoms_and_residues(tmp_path):
    path = _structure(tmp_path, "model.pdb", "pdb")
    kpi = CAtomCountPerformance()
    kpi.setFromPdbDataFile(_FileStub(path))
    assert int(kpi.nAtoms) == 3
    assert int(kpi.nResidues) == 2


def test_accepts_a_bare_path(tmp_path):
    """Some callers pass the object; str() of it is the path either way."""
    path = _structure(tmp_path, "model.pdb", "pdb")
    kpi = CAtomCountPerformance()
    kpi.setFromPdbDataFile(path)
    assert int(kpi.nAtoms) == 3


def test_mmcif_is_counted_the_same(tmp_path):
    """The format is read from the content, so both give the same answer."""
    pdb = CAtomCountPerformance()
    pdb.setFromPdbDataFile(_structure(tmp_path, "model.pdb", "pdb"))
    cif = CAtomCountPerformance()
    cif.setFromPdbDataFile(_structure(tmp_path, "model.cif", "mmcif"))
    assert (int(pdb.nAtoms), int(pdb.nResidues)) == (int(cif.nAtoms), int(cif.nResidues))


def test_a_missing_file_is_an_error_not_a_zero(tmp_path):
    """Silently reporting 0 atoms would be worse than failing."""
    kpi = CAtomCountPerformance()
    with pytest.raises(Exception):
        kpi.setFromPdbDataFile(str(tmp_path / "nope.pdb"))
