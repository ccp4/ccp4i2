"""MrParse hits from one entry matching several targets make a complex
template (docs/multi-component-mr.md, route A).

Built on the beta-lactamase / BLIP demo: two "cuts" (the demo models given
chain ids and a cell, as MrParse leaves an entry's chains) and a hits list
as homologs.json records them, with query ranges on the merged beta+BLIP
sequence.
"""
import json
from pathlib import Path

import gemmi

from ccp4i2.lib.utils.formats.mrparse_complexes import (
    describe, find_complexes, read_targets, target_of, write_complex,
)

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _seq(name):
    return "".join(l.strip() for l in (DEMO / f"{name}.seq").read_text().splitlines()
                   if not l.startswith(">"))


def _cut(path, name, chain_id):
    """A chain of an entry as MrParse cuts it: the entry's cell and frame."""
    structure = gemmi.read_structure(str(DEMO / f"{name}.pdb"))
    structure[0][0].name = chain_id
    structure.cell = gemmi.UnitCell(70.0, 80.0, 90.0, 90.0, 90.0, 90.0)
    structure.spacegroup_hm = "P 21 21 21"
    structure.setup_entities()
    structure.write_pdb(str(path))
    return path


def _hit(entry, chain, start, stop, identity, pdb_file):
    return {"name": f"{entry}_{chain}", "pdb_id": entry, "chain_id": chain,
            "query_start": start, "query_stop": stop, "seq_ident": identity,
            "pdb_file": f"homologs/{pdb_file}"}


BETA, BLIP = _seq("beta"), _seq("blip")
TARGETS = [("BETA", BETA), ("BLIP", BLIP)]
# On the merged sequence BETA (1-263) comes first, then BLIP (264-428)
HITS = [
    _hit("1xxx", "A", 3, 260, 0.95, "1xxx_A_3-260.pdb"),
    _hit("1xxx", "B", 266, 428, 0.80, "1xxx_B_2-165.pdb"),
    _hit("2yyy", "A", 1, 263, 1.00, "2yyy_A_1-263.pdb"),      # beta alone, another entry
    _hit("3zzz", "A", 5, 250, 0.60, "3zzz_A_5-250.pdb"),      # two chains of one entry,
    _hit("3zzz", "B", 7, 255, 0.60, "3zzz_B_7-255.pdb"),      # both beta: a homodimer, not a complex
]


def test_targets_are_read_as_mrparse_merges_them(tmp_path):
    fasta = tmp_path / "two.fasta"
    fasta.write_text(f">BETA some words\n{BETA[:30]}\n{BETA[30:]}\n>BLIP\n{BLIP}\n>BETA again\n{BETA}\n")
    assert read_targets(fasta) == TARGETS  # the repeated sequence is dropped


def test_a_hit_lies_on_the_target_holding_its_midpoint():
    assert target_of(HITS[0], TARGETS) == 0
    assert target_of(HITS[1], TARGETS) == 1
    assert target_of({"query_start": "x"}, TARGETS) is None


def test_only_an_entry_matching_different_targets_is_a_complex():
    found = find_complexes(HITS, TARGETS)
    assert [f["entry"] for f in found] == ["1xxx"]
    assert [(c["chain"], c["target"], c["identity"]) for c in found[0]["chains"]] == \
        [("A", "BETA", 0.95), ("B", "BLIP", 0.80)]
    assert describe(found[0]) == "1xxx chains A (BETA 95%), B (BLIP 80%)"


def test_the_template_holds_both_chains_in_the_entry_frame(tmp_path):
    cuts = tmp_path / "homologs"
    cuts.mkdir()
    _cut(cuts / "1xxx_A_3-260.pdb", "beta", "A")
    _cut(cuts / "1xxx_B_2-165.pdb", "blip", "B")
    entry = find_complexes(HITS, TARGETS)[0]
    out = write_complex(entry, cuts, tmp_path)
    assert out is not None and out.endswith("1xxx_AB_complex.pdb")
    structure = gemmi.read_structure(out)
    structure.setup_entities()
    assert [(c.name, len(c.get_polymer())) for c in structure[0]] == [("A", 263), ("B", 165)]
    assert structure.cell.a == 70.0 and structure.spacegroup_hm == "P 21 21 21"
    # the chains are the cuts' coordinates, unmoved
    beta = gemmi.read_structure(str(cuts / "1xxx_A_3-260.pdb"))[0][0][0][0].pos
    assert structure[0]["A"][0][0].pos.dist(beta) == 0.0


def test_a_missing_cut_makes_no_template(tmp_path):
    entry = find_complexes(HITS, TARGETS)[0]
    assert write_complex(entry, tmp_path / "nowhere", tmp_path) is None


def test_homologs_json_records_are_what_the_wrapper_passes(tmp_path):
    # the fields used are those of a real homologs.json (6p8e, CDK4/cyclin D1)
    record = {"name": "6p8e_B", "pdb_id": "6p8e", "chain_id": "B", "region_id": 1,
              "query_start": 18, "query_stop": 299, "seq_ident": 0.92,
              "pdb_file": "homologs/6p8e_B_19-297.pdb", "ellg": 1191.5}
    other = dict(record, name="6p8e_A", chain_id="A", query_start=322, query_stop=568,
                 seq_ident=1.0, pdb_file="homologs/6p8e_A_20-266.pdb")
    targets = [("CDK4", "M" * 303), ("CCND1", "M" * 295)]
    found = find_complexes([record, other], targets)
    assert describe(found[0]) == "6p8e chains B (CDK4 92%), A (CCND1 100%)"
    json.dumps(found)  # plain data, as program.xml and the report need
