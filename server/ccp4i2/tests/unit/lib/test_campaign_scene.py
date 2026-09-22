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


def test_placeholder_code_without_dictionary(tmp_path):
    """Without a dictionary, a placeholder LIG still registers."""
    coords = _write_pdb(tmp_path / "model.pdb", ligand_code="LIG")
    assert campaign_scene.detect_ligands(coords, dict_path=None) == ["LIG"]


def test_real_ligand_code_detected_without_dictionary(tmp_path):
    """A real soaked-fragment code (e.g. NUT) is detected even with no dict.

    Fragment campaigns refined with servalcat often carry no restraint
    dictionary, so detection must not depend on one.
    """
    coords = _write_pdb(tmp_path / "model.pdb", ligand_code="NUT")
    assert campaign_scene.detect_ligands(coords, dict_path=None) == ["NUT"]


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


# --------------------------------------------------------------------------
# Site geometry: the pocket and the nearest-fragment diagnostic
# --------------------------------------------------------------------------

def _residue(name, seqnum, atoms, icode=" "):
    """A residue with named atoms at given positions."""
    res = gemmi.Residue()
    res.name = name
    res.seqid = gemmi.SeqId(seqnum, icode)
    for atom_name, (x, y, z) in atoms:
        atom = gemmi.Atom()
        atom.name = atom_name
        atom.element = gemmi.Element("O" if name == "HOH" else "C")
        atom.pos = gemmi.Position(x, y, z)
        res.add_atom(atom)
    return res


def _structure(chains):
    """``{chain_name: [residues]}`` -> a one-model structure."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(100, 100, 100, 90, 90, 90)
    model = gemmi.Model(1)
    for chain_name, residues in chains.items():
        chain = gemmi.Chain(chain_name)
        for res in residues:
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    return st


ORIGIN = (0.0, 0.0, 0.0)


def test_pocket_residues_within_radius_as_cids():
    """A residue is in the pocket if ANY atom is in range, and appears once.

    Residue 2 has its CA out of range and one side-chain atom inside it, so
    it belongs; residue 3 is out of range entirely.
    """
    st = _structure({
        "A": [
            _residue("ALA", 1, [("CA", (2, 0, 0)), ("CB", (3, 0, 0))]),
            _residue("LYS", 2, [("CA", (12, 0, 0)), ("NZ", (7, 0, 0))]),
            _residue("GLY", 3, [("CA", (20, 0, 0))]),
        ]
    })
    assert campaign_scene.pocket_residue_cids(st, ORIGIN, radius=8.0) == [
        "//A/1", "//A/2",
    ]


def test_pocket_is_in_chain_and_sequence_order_not_lexical():
    """``//A/100`` sorts after ``//A/45`` numerically; a string sort would not."""
    st = _structure({
        "B": [_residue("ALA", 7, [("CA", (1, 0, 0))])],
        "A": [
            _residue("ALA", 100, [("CA", (0, 1, 0))]),
            _residue("ALA", 45, [("CA", (0, 0, 1))]),
        ],
    })
    assert campaign_scene.pocket_residue_cids(st, ORIGIN) == [
        "//A/45", "//A/100", "//B/7",
    ]


def test_pocket_insertion_code_follows_a_dot():
    """mmdb reads ``45.A``; ``45A`` would select nothing."""
    st = _structure({
        "A": [_residue("ALA", 45, [("CA", (1, 0, 0))], icode="A")]
    })
    assert campaign_scene.pocket_residue_cids(st, ORIGIN) == ["//A/45.A"]


def test_pocket_in_empty_solvent_is_empty():
    """The list must be empty, so the builder can emit no representation
    at all; an empty selection would draw the whole molecule."""
    st = _structure({"A": [_residue("ALA", 1, [("CA", (2, 0, 0))])]})
    assert campaign_scene.pocket_residue_cids(st, (50, 50, 50)) == []


def test_pocket_leaves_out_water_and_fragments():
    """The pocket is what the ligands bind to: not the solvent, and not a
    fragment (which, when the exemplar is a hit, is drawn as a hit)."""
    st = _structure({
        "A": [
            _residue("ALA", 1, [("CA", (2, 0, 0))]),
            _residue("HOH", 201, [("O", (1, 0, 0))]),
            _residue("DRG", 301, [("C1", (0, 1, 0))]),
            _residue("SO4", 401, [("S", (0, 0, 1))]),
        ]
    })
    # A sulphate is context, not a fragment, so it stays.
    assert campaign_scene.pocket_residue_cids(st, ORIGIN) == ["//A/1", "//A/401"]


def test_nearest_fragment_distance_is_to_the_closest_centroid():
    st = _structure({
        "A": [
            _residue("ALA", 1, [("CA", (0.5, 0, 0))]),
            _residue("DRG", 301, [("C1", (9, 0, 0)), ("C2", (11, 0, 0))]),
            _residue("LIG", 302, [("C1", (0, 30, 0))]),
        ]
    })
    assert campaign_scene.nearest_fragment_distance(st, ORIGIN) == pytest.approx(10.0)


def test_nearest_fragment_distance_is_none_without_a_fragment():
    """None, not inf: a stats payload with ``inf`` in it is not JSON."""
    st = _structure({
        "A": [
            _residue("ALA", 1, [("CA", (0.5, 0, 0))]),
            _residue("HOH", 201, [("O", (1, 0, 0))]),
        ]
    })
    assert campaign_scene.nearest_fragment_distance(st, ORIGIN) is None


def test_nearest_fragment_distance_measures_after_the_transform():
    """Measured in the frame the scene draws: a fit that brings the
    fragment onto the site makes the distance small."""
    st = _structure({
        "A": [_residue("DRG", 301, [("C1", (10, 0, 0))])]
    })
    shift = gemmi.Transform(gemmi.Mat33(), gemmi.Vec3(-10, 0, 0))
    assert campaign_scene.nearest_fragment_distance(st, ORIGIN) == pytest.approx(10.0)
    assert campaign_scene.nearest_fragment_distance(st, ORIGIN, shift) == pytest.approx(0.0)


# --------------------------------------------------------------------------
# site_position: the sign that cost a whole implementation round
# --------------------------------------------------------------------------

class _FakeSite:
    def __init__(self, x, y, z):
        self.origin_x, self.origin_y, self.origin_z = x, y, z

    @property
    def origin(self):
        return [self.origin_x, self.origin_y, self.origin_z]


def test_site_position_negates_the_stored_origin():
    """A CampaignSite stores Moorhen's view origin, not a position.

    Moorhen's origin is the negation of the point at screen centre, and the
    site save/restore path stores and restores it raw, so the two negations
    cancel and nothing looks wrong until something treats the stored value as
    a coordinate. The first thing that did -- the site scene's pocket
    selection -- looked 41 A from the protein and found nothing.
    """
    site = _FakeSite(24.303, 8.865, -0.941)
    assert campaign_scene.site_position(site) == (-24.303, -8.865, 0.941)


def test_site_position_is_what_finds_the_pocket(tmp_path):
    """The negated origin selects residues; the stored one selects nothing.

    Pins the direction, not just the arithmetic: a sign flip that still
    negated *something* would pass the test above and fail this one.
    """
    path = tmp_path / "pocket.pdb"
    _write_pdb(path, ligand_code=None)
    st = gemmi.read_structure(str(path))

    centre = campaign_scene.pocket_residue_cids(st, (0.0, 0.0, 0.0), radius=20.0)
    assert centre, "the fixture should have residues near the real origin"

    site = _FakeSite(0.0, 0.0, 0.0)
    assert campaign_scene.pocket_residue_cids(
        st, campaign_scene.site_position(site), radius=20.0
    ) == centre

    # Off-centre: stored (30,0,0) means the screen was centred on (-30,0,0).
    away = _FakeSite(300.0, 0.0, 0.0)
    assert campaign_scene.pocket_residue_cids(
        st, campaign_scene.site_position(away), radius=20.0
    ) == []
