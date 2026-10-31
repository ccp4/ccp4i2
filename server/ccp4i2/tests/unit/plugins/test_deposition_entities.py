"""The deposited model carries one entity per crystallised molecule, with its full sequence.

Refinement programs write one entity per chain and no sequence. The PDBe
adding_stats_to_mmcif package, asked to add the sequences, appends
_entity_poly rows one entity at a time and re-reads the loop between appends;
on the Ahir3252 hexamer (six servalcat entities) that produced a 23,000-line
_entity_poly loop with the sequence in the wrong column, which OneDep could
not read. So the deposition task now builds the entities itself from the AU
contents, before the package runs, and restores its table afterwards.

These tests use a two-chain toy: chain A and chain B are the same molecule,
each missing residues at the ends and one in the middle, as a real model does.
"""
import gemmi
import pytest

from ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.deposition_entities import (
    entity_summary,
    rebuild_entities,
    rewrite_entity_poly,
)

# The expressed construct: the model below starts at residue 3 and lacks residue 6.
CONSTRUCT = "MGSAKLVADER"
MODELLED = {3: "SER", 4: "ALA", 5: "LYS", 7: "VAL", 8: "ALA", 9: "ASP"}


def _refinement_style_mmcif(path, chain_b=MODELLED):
    """What servalcat/refmac hand over: per-chain entities, '?' sequences, a _refine table, and a
    _pdbx_database_status copied from whatever PDB file the model came from."""
    st = gemmi.Structure()
    st.name = "refined"
    st.cell = gemmi.UnitCell(50, 60, 70, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    for chain_id, residues in (("A", MODELLED), ("B", chain_b)):
        chain = gemmi.Chain(chain_id)
        for num, name in residues.items():
            res = gemmi.Residue()
            res.name, res.seqid = name, gemmi.SeqId(num, " ")
            for i, atom_name in enumerate(("N", "CA", "C", "O")):
                atom = gemmi.Atom()
                atom.name, atom.element = atom_name, gemmi.Element("C" if atom_name != "N" else "N")
                atom.pos = gemmi.Position(num * 3.8 + i, 10 if chain_id == "A" else 30, 0)
                atom.occ, atom.b_iso = 1.0, 30.0
                res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    doc = st.make_mmcif_document()
    block = doc[0]
    # one entity per chain, no sequence: the refinement programs' habit
    block.set_mmcif_category("_entity_poly.", {
        "entity_id": ["A", "B"], "type": ["polypeptide(L)"] * 2,
        "pdbx_strand_id": ["A", "B"], "pdbx_seq_one_letter_code": ["?", "?"]})
    block.set_mmcif_category("_refine.", {"entry_id": ["XXXX"], "ls_d_res_high": ["2.1"]})
    block.set_mmcif_category("_pdbx_database_status.", {
        "entry_id": ["XXXX"], "recvd_initial_deposition_date": ["2015-12-22"]})
    doc.write_file(str(path))


@pytest.fixture
def refined(tmp_path):
    path = tmp_path / "refined.cif"
    _refinement_style_mmcif(path)
    return path


def test_one_entity_with_the_full_sequence(refined, tmp_path):
    out = tmp_path / "rebuilt.cif"
    rebuild_entities(refined, out, [("construct", CONSTRUCT, "PROTEIN")])
    n, strands, lengths = entity_summary(out)
    assert n == 1
    assert sorted(strands["1"]) == ["A", "B"]
    assert lengths["1"] == len(CONSTRUCT)
    block = gemmi.cif.read(str(out))[0]
    seq = block.find("_entity_poly_seq.", ["entity_id", "num", "mon_id"])
    assert len(seq) == len(CONSTRUCT)
    assert [row.str(2) for row in seq][:3] == ["MET", "GLY", "SER"]


def test_label_seq_id_follows_the_construct(refined, tmp_path):
    out = tmp_path / "rebuilt.cif"
    rebuild_entities(refined, out, [("construct", CONSTRUCT, "PROTEIN")])
    block = gemmi.cif.read(str(out))[0]
    seen = {}
    tags = ["auth_asym_id", "auth_seq_id", "label_seq_id", "label_entity_id"]
    for row in block.find("_atom_site.", tags):
        seen[(row.str(0), int(row.str(1)))] = (int(row.str(2)), row.str(3))
    # residue 7 of the model is residue 7 of the construct, in entity 1, in both chains
    assert seen[("A", 7)] == (7, "1")
    assert seen[("B", 7)] == (7, "1")
    assert seen[("A", 3)] == (3, "1")


def test_refinement_categories_are_kept_and_the_borrowed_header_is_not(refined, tmp_path):
    out = tmp_path / "rebuilt.cif"
    rebuild_entities(refined, out, [("construct", CONSTRUCT, "PROTEIN")])
    block = gemmi.cif.read(str(out))[0]
    assert block.find_value("_refine.ls_d_res_high") == "2.1"
    assert not block.find_mmcif_category("_pdbx_database_status.")


def test_two_sequences_each_get_their_chains(tmp_path):
    """A complex: chain B is a different molecule from chain A."""
    path = tmp_path / "refined.cif"
    # chain B is a different protein (Lys->Trp, Asp->Phe)
    chain_b = {num: {"LYS": "TRP", "ASP": "PHE"}.get(name, name) for num, name in MODELLED.items()}
    _refinement_style_mmcif(path, chain_b)
    out = tmp_path / "rebuilt.cif"
    sequences = [("alpha", CONSTRUCT, "PROTEIN"), ("beta", "MGSAWLVAFER", "PROTEIN")]
    rebuild_entities(path, out, sequences)
    n, strands, lengths = entity_summary(out)
    assert n == 2
    assert strands == {"1": ["A"], "2": ["B"]}


def test_rewrite_entity_poly_undoes_the_package_append(refined, tmp_path):
    out = tmp_path / "rebuilt.cif"
    st = rebuild_entities(refined, out, [("construct", CONSTRUCT, "PROTEIN")])
    # what adding_stats_to_mmcif leaves behind: the row appended, type column blanked
    doc = gemmi.cif.read(str(out))
    doc[0].set_mmcif_category("_entity_poly.", {
        "entity_id": ["1", "1", "1"], "type": ["polypeptide(L)", "polypeptide(L)", ""],
        "pdbx_strand_id": ["", "", ""], "pdbx_seq_one_letter_code": [CONSTRUCT] * 3})
    doc.write_file(str(out))
    rewrite_entity_poly(out, st)
    n, strands, lengths = entity_summary(out)
    assert (n, sorted(strands["1"]), lengths["1"]) == (1, ["A", "B"], len(CONSTRUCT))


def test_no_sequence_is_refused(refined, tmp_path):
    with pytest.raises(ValueError):
        rebuild_entities(refined, tmp_path / "out.cif", [])
