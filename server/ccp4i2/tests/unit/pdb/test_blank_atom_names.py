"""
Regression tests for atoms written without a name.

Moorhen's "delete residue, add it back as ALA" creates the backbone N/CA/C/O
but leaves the carbonyl oxygen's name blank while type_symbol is correctly O.
An unnamed atom matches no monomer-library definition, so refinement -- and
CCP4i2's own checkMonomeCoverage pre-flight -- reported the whole residue as
having atoms without restraints, with the message mangled by the empty name
("(replace  with N)").  Two real files from one session carried 9 and 5 such
atoms; the coordinates themselves were sound.

CPdbData.loadFile must therefore recover the name from the element where that
is unambiguous, and must not invent names it cannot know.
"""

import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi (CCP4 python)")

from ccp4i2.core.CCP4ModelData import CPdbData, _repair_blank_atom_names


def _residue(name, seq, atoms):
    res = gemmi.Residue()
    res.name = name
    res.seqid = gemmi.SeqId(str(seq))
    for nm, el in atoms:
        a = gemmi.Atom()
        a.name = nm
        a.element = gemmi.Element(el)
        a.pos = gemmi.Position(1, 2, 3)
        a.occ = 1.0
        a.b_iso = 20.0
        res.add_atom(a)
    return res


def _structure(residues):
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")
    for res in residues:
        chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    return st


def test_blank_carbonyl_oxygen_is_named_from_element():
    """The observed Moorhen defect: O written with no name."""
    st = _structure([
        _residue("GLY", 15, [("", "O"), ("N", "N"), ("CA", "C"), ("C", "C")]),
    ])
    repaired = _repair_blank_atom_names(st)

    assert len(repaired) == 1
    names = sorted(a.name for a in st[0]["A"][0])
    assert names == ["C", "CA", "N", "O"]


def test_whitespace_only_name_is_treated_as_blank():
    """A name of pure padding is as useless as an empty one."""
    st = _structure([
        _residue("ALA", 16, [("   ", "O"), ("N", "N"), ("CA", "C"), ("C", "C")]),
    ])
    assert len(_repair_blank_atom_names(st)) == 1
    assert "O" in {a.name for a in st[0]["A"][0]}


def test_blank_name_not_overwritten_when_residue_already_has_that_atom():
    """Never create a duplicate: if O is present, leave the blank alone."""
    st = _structure([
        _residue("GLY", 17, [("O", "O"), ("", "O"), ("N", "N"), ("CA", "C")]),
    ])
    assert _repair_blank_atom_names(st) == []
    assert sum(1 for a in st[0]["A"][0] if not a.name.strip()) == 1


def test_ambiguous_elements_are_left_alone():
    """C and N appear many times in a residue; guessing would be wrong."""
    st = _structure([
        _residue("LEU", 22, [("", "C"), ("", "N"), ("N", "N"), ("CA", "C")]),
    ])
    assert _repair_blank_atom_names(st) == []
    assert sum(1 for a in st[0]["A"][0] if not a.name.strip()) == 2


def test_wellformed_structure_is_untouched():
    st = _structure([
        _residue("ALA", 1, [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]),
    ])
    assert _repair_blank_atom_names(st) == []


def test_loadfile_repairs_blank_names(tmp_path):
    """End to end: the repair happens on load, so nothing downstream sees it."""
    src = tmp_path / "blank.cif"
    st = _structure([
        _residue("GLY", 15, [("", "O"), ("N", "N"), ("CA", "C"), ("C", "C")]),
        _residue("ALA", 16, [("", "O"), ("N", "N"), ("CA", "C"), ("C", "C"), ("CB", "C")]),
    ])
    st.make_mmcif_document().write_file(str(src))

    # Sanity: the fixture really does round-trip with blank names.
    reread = gemmi.read_structure(str(src))
    assert sum(1 for ch in reread[0] for r in ch for a in r if not a.name.strip()) == 2

    pdb = CPdbData()
    err = pdb.loadFile(str(src))
    assert err.count() == 0, f"loadFile errored: {err}"

    loaded = pdb._gemmi_structure
    assert sum(1 for ch in loaded[0] for r in ch for a in r if not a.name.strip()) == 0
    for res in loaded[0]["A"]:
        assert "O" in {a.name for a in res}


def test_repair_preserves_coordinates_occupancy_and_b():
    """Only the name may change -- the atom is otherwise sound."""
    st = _structure([
        _residue("GLY", 15, [("", "O"), ("N", "N"), ("CA", "C"), ("C", "C")]),
    ])
    atom = st[0]["A"][0][0]
    before = (atom.pos.x, atom.pos.y, atom.pos.z, atom.occ, atom.b_iso, atom.element.name)

    _repair_blank_atom_names(st)

    atom = [a for a in st[0]["A"][0] if a.name == "O"][0]
    assert (atom.pos.x, atom.pos.y, atom.pos.z, atom.occ, atom.b_iso,
            atom.element.name) == before


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
