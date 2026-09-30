"""The monomer payload carries what a 2D depiction needs, in a fixed order.

A picker built on this renders its own SVG so that every element carries the
dictionary atom name ("NZ"), which is the currency the AceDRG tasks speak.
That only works if the atom list order is a contract: index i means atoms[i],
so a molecule built from this never has to map depiction atom indices back to
names.
"""
import pytest

from ccp4i2.lib.utils.formats.cif_ligand import extract_monomer_atoms_bonds

gemmi = pytest.importorskip("gemmi", reason="the extractor is gemmi-based")

# A two-heavy-atom monomer with a hydrogen (dropped) and a formal charge.
CIF = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
XYZ XYZ 'test' non-polymer 4 3

data_comp_XYZ
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
XYZ N1  N  1.000  0.000 0.000 0.000
XYZ C2  C  0.000  1.500 0.000 0.000
XYZ O3  O  -1.000 2.200 1.000 0.000
XYZ H1  H  0.000 -0.500 0.800 0.000
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
XYZ N1 C2 single 1.500 0.020
XYZ C2 O3 double 1.230 0.020
XYZ N1 H1 single 1.000 0.020
"""


@pytest.fixture
def monomer(tmp_path):
    path = tmp_path / "XYZ.cif"
    path.write_text(CIF)
    return extract_monomer_atoms_bonds(str(path))


def test_atom_details_track_atoms_exactly(monomer):
    # The order contract: same atoms, same order, so index i means atoms[i].
    assert monomer["atoms"] == [d["name"] for d in monomer["atom_details"]]


def test_hydrogens_are_dropped_from_both_lists(monomer):
    assert "H1" not in monomer["atoms"]
    assert "H1" not in [d["name"] for d in monomer["atom_details"]]
    # ... and so is the bond that used one.
    assert all("H1" not in (b["atom1"], b["atom2"]) for b in monomer["bonds"])


def test_elements_are_carried(monomer):
    assert {d["name"]: d["element"] for d in monomer["atom_details"]} == {
        "N1": "N", "C2": "C", "O3": "O"}


def test_formal_charge_is_a_whole_number(monomer):
    # Dictionaries write these as floats ("1.000"); a molfile wants an int.
    charges = {d["name"]: d["charge"] for d in monomer["atom_details"]}
    assert charges == {"N1": 1, "C2": 0, "O3": -1}
    assert all(isinstance(value, int) for value in charges.values())


def test_hydrogens_are_counted_on_the_atom_that_carries_them(monomer):
    # Not drawn, but a valence count over an edited monomer needs them.
    assert {d["name"]: d["hydrogens"] for d in monomer["atom_details"]} == {
        "N1": 1, "C2": 0, "O3": 0}


def test_bond_orders_survive(monomer):
    assert {(b["atom1"], b["atom2"]): b["type"] for b in monomer["bonds"]} == {
        ("N1", "C2"): "single", ("C2", "O3"): "double"}


def test_an_unreadable_file_yields_the_same_shape(tmp_path):
    junk = tmp_path / "junk.cif"
    junk.write_text("not a cif at all\n")
    result = extract_monomer_atoms_bonds(str(junk))
    # Callers index all three keys; none may be missing on the failure path.
    assert result == {"atoms": [], "bonds": [], "atom_details": []}
