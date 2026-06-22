"""
Unit tests for fragment-campaign hit detection (``lib.campaign_scene``).

These exercise the pure, gemmi-only detection helpers — no Django models,
no CCP4 binaries — using tiny coordinate + dictionary files built in a
tmp dir. The scene-assembly half (``build_summary_scene``) touches the DB
and is covered by the API e2e suite instead.
"""

import gemmi
import pytest

from ccp4i2.lib import campaign_scene


# --------------------------------------------------------------------------
# Fixtures: build minimal coord + dict files
# --------------------------------------------------------------------------

def _write_pdb(path, ligand_code=None):
    """A one-residue alanine 'protein', optionally plus a HET ligand."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")

    res = gemmi.Residue()
    res.name = "ALA"
    res.seqid = gemmi.SeqId("1")
    for nm, el in [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]:
        a = gemmi.Atom()
        a.name = nm
        a.element = gemmi.Element(el)
        a.pos = gemmi.Position(1, 2, 3)
        a.occ = 1.0
        a.b_iso = 20.0
        res.add_atom(a)
    chain.add_residue(res)

    if ligand_code:
        lig = gemmi.Residue()
        lig.name = ligand_code
        lig.seqid = gemmi.SeqId("2")
        lig.het_flag = "H"
        a = gemmi.Atom()
        a.name = "C1"
        a.element = gemmi.Element("C")
        a.pos = gemmi.Position(5, 5, 5)
        a.occ = 1.0
        a.b_iso = 30.0
        lig.add_atom(a)
        chain.add_residue(lig)

    model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    path.write_text(st.make_pdb_string())
    return path


def _write_dict(path, *codes):
    """A minimal refmac-style restraint CIF declaring one or more monomers.

    Passing several codes simulates a merged LIBOUT that carries standard
    monomers (ALA, GLY, ...) alongside the actual ligand.
    """
    list_rows = "\n".join(
        f"{c} {c} 'monomer' non-polymer 3 3 ." for c in codes
    )
    comp_blocks = "\n".join(
        f"""data_comp_{c}
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
{c} C1 C
{c} C2 C
{c} O1 O"""
        for c in codes
    )
    text = f"""data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
{list_rows}
{comp_blocks}
"""
    path.write_text(text)
    return path


# --------------------------------------------------------------------------
# dictionary_comp_ids / coordinate_residue_names
# --------------------------------------------------------------------------

def test_dictionary_comp_ids_reads_all_sources(tmp_path):
    d = _write_dict(tmp_path / "DRG.cif", "DRG")
    ids = campaign_scene.dictionary_comp_ids(d)
    assert "DRG" in ids
    # The catch-all list block must not leak in as a comp id.
    assert "LIST" not in ids


def test_coordinate_residue_names(tmp_path):
    p = _write_pdb(tmp_path / "model.pdb", ligand_code="DRG")
    names = campaign_scene.coordinate_residue_names(p)
    assert "ALA" in names
    assert "DRG" in names


# --------------------------------------------------------------------------
# detect_ligands — the hit decision
# --------------------------------------------------------------------------

def test_dictionary_driven_hit(tmp_path):
    """Ligand named for its dict comp_id is detected even with a non-LIG code."""
    coords = _write_pdb(tmp_path / "model.pdb", ligand_code="ABC")
    dictf = _write_dict(tmp_path / "ABC.cif", "ABC")
    assert campaign_scene.detect_ligands(coords, dictf) == ["ABC"]


def test_dictionary_present_but_ligand_absent_is_not_a_hit(tmp_path):
    """An apo dataset (dict exists, but no such residue in coords) is no hit."""
    coords = _write_pdb(tmp_path / "apo.pdb", ligand_code=None)
    dictf = _write_dict(tmp_path / "DRG.cif", "DRG")
    assert campaign_scene.detect_ligands(coords, dictf) == []


def test_fallback_codes_when_no_dictionary(tmp_path):
    """Without a dictionary, LIG/DRG/UNL placeholder codes still register."""
    coords = _write_pdb(tmp_path / "model.pdb", ligand_code="LIG")
    assert campaign_scene.detect_ligands(coords, dict_path=None) == ["LIG"]


def test_fallback_ignores_unknown_code_without_dictionary(tmp_path):
    """A real 3-letter code is invisible to the no-dictionary fallback."""
    coords = _write_pdb(tmp_path / "model.pdb", ligand_code="ABC")
    assert campaign_scene.detect_ligands(coords, dict_path=None) == []


def test_merged_dict_with_standard_monomers_does_not_false_positive(tmp_path):
    """A LIBOUT carrying ALA/GLY alongside the ligand flags only the ligand."""
    coords = _write_pdb(tmp_path / "model.pdb", ligand_code="ABC")
    # Dict declares the fragment plus standard amino acids (as merged LIBOUTs do).
    dictf = _write_dict(tmp_path / "merged.cif", "ABC", "ALA", "GLY")
    assert campaign_scene.detect_ligands(coords, dictf) == ["ABC"]


def test_apo_with_merged_dict_is_not_a_hit(tmp_path):
    """An apo model whose only dict comps are standard residues is no hit."""
    coords = _write_pdb(tmp_path / "apo.pdb", ligand_code=None)
    dictf = _write_dict(tmp_path / "merged.cif", "DRG", "ALA", "GLY")
    assert campaign_scene.detect_ligands(coords, dictf) == []


def test_crystallisation_additive_is_not_a_hit(tmp_path):
    """Glycerol bound (and in the dict) must not register as a fragment."""
    coords = _write_pdb(tmp_path / "model.pdb", ligand_code="GOL")
    dictf = _write_dict(tmp_path / "GOL.cif", "GOL")
    assert campaign_scene.detect_ligands(coords, dictf) == []


def test_unreadable_coordinates_are_not_a_hit(tmp_path):
    bad = tmp_path / "nope.pdb"
    bad.write_text("this is not a coordinate file\n")
    # Should swallow the gemmi error and report no hit, not raise.
    assert campaign_scene.detect_ligands(bad, dict_path=None) == []
