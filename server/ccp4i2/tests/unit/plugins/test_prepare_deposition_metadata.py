"""The descriptive metadata a deposition carries is written into the model
mmCIF, every category once.

OneDep reads title, authors, source organisms, sequence references, data
collection and assembly from the model file; a category written in two
places is rejected as "Bad model file" (learned at the cost of eight uploads
on 2026-11-01). These tests drive ccp4i2.lib.deposition.metadata on a toy
two-chain model through the four kinds of molecule the task must describe:
a natural protein with an engineered variant, a de novo design, a chemically
synthesised peptide, and a chimera of designed and natural segments.
"""
import os

import gemmi
import pytest

from ccp4i2.lib.deposition import metadata, mmcif_edit, statistics
from ccp4i2.lib.deposition.entities import rebuild_entities_in
from ccp4i2.lib.deposition.struct_ref import (
    DIF_CONFLICT,
    DIF_ENGINEERED,
    DIF_INSERTION,
    DIF_TAG,
    Reference,
    parse_fasta,
)

DICTIONARY = os.path.join(os.environ.get("CCP4", ""), "share", "libcifpp", "mmcif_pdbx.dic")

# A natural protein (the "UniProt" entry) and the construct made from it: an
# N-terminal tag, two point changes (D8H intended, K2Q a cloning artefact),
# a three-residue insertion, and a C-terminal cloning remnant.
REFERENCE = "MKIYYEHDLDGFILGKKVAVLG"
TAG = "MGHHHHHH"
CORE = "MQIYYEHHLDGFILGGGGKKVAVLG"   # K2Q, D8H, GGG inserted after G15
TAIL = "RS"
CONSTRUCT = TAG + CORE + TAIL


def _toy_model(path, construct, chains=("A", "B"), start=3, skip=()):
    """A refinement-style mmCIF: each chain models ``construct`` from
    residue ``start`` with the residues in ``skip`` missing."""
    st = gemmi.Structure()
    st.name = "refined"
    st.cell = gemmi.UnitCell(50, 60, 70, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    for ci, chain_id in enumerate(chains):
        chain = gemmi.Chain(chain_id)
        for num, one in enumerate(construct, start=1):
            if num < start or num in skip:
                continue
            res = gemmi.Residue()
            res.name = gemmi.find_tabulated_residue(
                gemmi.expand_one_letter(one, gemmi.ResidueKind.AA)).name
            res.seqid = gemmi.SeqId(num, " ")
            for i, atom_name in enumerate(("N", "CA", "C", "O")):
                atom = gemmi.Atom()
                atom.name, atom.element = atom_name, gemmi.Element("N" if atom_name == "N" else "C")
                atom.pos = gemmi.Position(num * 3.8 + i, 10 + 20 * ci, 0)
                atom.occ, atom.b_iso = 1.0, 30.0
                res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    doc = st.make_mmcif_document()
    block = doc[0]
    block.set_mmcif_category("_refine.", {"entry_id": ["XXXX"], "ls_d_res_high": ["2.1"],
                                          "ls_d_res_low": ["40.0"], "ls_R_factor_R_work": ["0.19"]})
    block.set_mmcif_category("_software.", {"pdbx_ordinal": ["1"], "name": ["refmac"],
                                            "version": ["5.8"], "classification": ["refinement"]})
    block.set_mmcif_category("_diffrn_radiation_wavelength.", {"id": ["1"], "wavelength": ["0.9795"]})
    doc.write_file(str(path))
    return path


def _prepared(tmp_path, construct=CONSTRUCT, **kw):
    doc = gemmi.cif.read(str(_toy_model(tmp_path / "model.cif", construct, **kw)))
    st = rebuild_entities_in(doc, [("construct", construct, "PROTEIN")])
    return doc, doc[0], st


def _fetch(accession):
    return "TOY_PROT", REFERENCE


BASE_SPEC = {
    "title": "A toy for the deposition tests",
    "keywords": "OXIDOREDUCTASE",
    "authors": [{"name": "Lewis, C."}, {"name": "Noble, M.E.M.", "orcid": "0000-0002-3595-9807"}],
    "contact": {"first_name": "Martin", "last_name": "Noble", "email": "m@example.org",
                "country": "United Kingdom"},
    "related": [{"db_name": "PDB", "db_id": "8pve", "details": "same protein"}],
    "experiment": {"site": "Diamond", "beamline": "I02", "detector": "PIXEL",
                   "detector_type": "DECTRIS PILATUS 6M-F", "collection_date": "2015-12-11",
                   "temperature": 100},
    "phasing": {"method": "MR", "starting_model": "2JAB"},
}


def _rows(block, category, tags):
    return [tuple(r.str(i) or "?" for i in range(len(tags))) for r in block.find(category + ".", tags)]


def test_every_category_is_written_once(tmp_path):
    doc, block, st = _prepared(tmp_path)
    spec = dict(BASE_SPEC, entities=[{"entity_id": "1", "description": "toy protein", "segments": [
        {"beg": 1, "end": 8, "kind": "tag"},
        {"beg": 9, "end": 33, "kind": "engineered", "organism": "Thermus thermophilus",
         "taxid": 274, "uniprot": "P00000", "engineered": "D8H"},
        {"beg": 34, "end": 35, "kind": "tag"}]}])
    warnings = metadata.write_metadata(block, st, spec, fetch=_fetch)
    assert warnings == []
    assert mmcif_edit.split_categories(block) == []
    # The two _refine items the writer adds join the existing _refine pairs.
    assert gemmi.cif.as_string(block.find_value("_refine.pdbx_starting_model")) == "2JAB"
    assert gemmi.cif.as_string(block.find_value("_refine.ls_R_factor_R_work")) == "0.19"
    # Written and re-read, still once each.
    out = tmp_path / "out.cif"
    doc.write_file(str(out), gemmi.cif.Style.Pdbx)
    assert mmcif_edit.split_categories(gemmi.cif.read(str(out))[0]) == []


@pytest.mark.skipif(not os.path.exists(DICTIONARY), reason="no PDBx dictionary in $CCP4")
def test_written_categories_conform_to_the_dictionary(tmp_path):
    doc, block, st = _prepared(tmp_path)
    spec = dict(BASE_SPEC, entities=[{"entity_id": "1", "description": "toy protein", "segments": [
        {"beg": 1, "end": 8, "kind": "tag"},
        {"beg": 9, "end": 33, "kind": "engineered", "organism": "Thermus thermophilus",
         "taxid": 274, "uniprot": "P00000", "engineered": "D8H"},
        {"beg": 34, "end": 35, "kind": "tag"}]}],
        assembly={"mode": "author", "oligomeric_count": 2, "operators": "-x,y+1/2,-z"})
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    messages = mmcif_edit.check_categories(block, DICTIONARY, metadata.WRITTEN_CATEGORIES)
    assert messages == [], messages


def test_engineered_segment_records_every_difference(tmp_path):
    doc, block, st = _prepared(tmp_path)
    spec = dict(BASE_SPEC, entities=[{"entity_id": "1", "description": "toy protein", "segments": [
        {"beg": 1, "end": 8, "kind": "tag"},
        {"beg": 9, "end": 33, "kind": "engineered", "organism": "Thermus thermophilus",
         "taxid": 274, "uniprot": "P00000", "engineered": "D8H"},
        {"beg": 34, "end": 35, "kind": "tag"}]}])
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    refs = _rows(block, "_struct_ref", ["id", "db_name", "db_code", "pdbx_db_accession", "entity_id"])
    assert refs == [("1", "UNP", "TOY_PROT", "P00000", "1")]
    seqs = _rows(block, "_struct_ref_seq", ["pdbx_strand_id", "seq_align_beg", "seq_align_end",
                                             "db_align_beg", "db_align_end"])
    # Two gapless stretches per chain: before and after the insertion (the
    # four glycines are ambiguous; the aligner takes the last as the gene's).
    assert sorted(seqs)[:2] == [("A", "26", "33", "15", "22"), ("A", "9", "22", "1", "14")]
    difs = {(r[0], r[1], r[2]): r[3] for r in _rows(
        block, "_struct_ref_seq_dif", ["pdbx_pdb_strand_id", "seq_num", "mon_id", "details"])}
    assert difs[("A", "10", "GLN")] == DIF_CONFLICT          # K2Q, not declared
    assert difs[("A", "16", "HIS")] == DIF_ENGINEERED        # D8H, declared
    assert difs[("A", "24", "GLY")] == DIF_INSERTION
    assert difs[("A", "1", "MET")] == DIF_TAG and difs[("A", "35", "SER")] == DIF_TAG
    gen = _rows(block, "_entity_src_gen", ["pdbx_src_id", "pdbx_seq_type", "pdbx_beg_seq_num",
                                           "pdbx_end_seq_num", "pdbx_gene_src_scientific_name",
                                           "pdbx_gene_src_ncbi_taxonomy_id"])
    assert gen[0] == ("1", "N-terminal tag", "1", "8", "synthetic construct", "32630")
    assert gen[1][1:] == ("Biological sequence", "9", "33", "Thermus thermophilus", "274")
    assert gen[2][1] == "C-terminal tag"
    entity = _rows(block, "_entity", ["id", "type", "src_method", "pdbx_description"])
    assert entity == [("1", "polymer", "man", "toy protein")]


def test_a_de_novo_design_references_the_entry_itself(tmp_path):
    design = "MSEEALKRLAEEAIKLAKEGNPEEARKLLEEAVRLAEEALK"
    doc, block, st = _prepared(tmp_path, construct=design)
    spec = dict(BASE_SPEC, keywords=None,
                design_details="RFdiffusion backbone, ProteinMPNN sequence, AF2 pLDDT 92",
                software=[{"name": "RFdiffusion", "version": "1.1", "classification": "model building"},
                          {"name": "ProteinMPNN", "classification": "model building"}],
                entities=[{"entity_id": "1", "description": "designed helical bundle",
                           "segments": [{"beg": 1, "end": len(design), "kind": "designed"}]}])
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    assert mmcif_edit.split_categories(block) == []
    assert gemmi.cif.as_string(block.find_value("_struct_keywords.pdbx_keywords")) == "DE NOVO PROTEIN"
    refs = _rows(block, "_struct_ref", ["db_name", "db_code", "pdbx_db_accession",
                                        "pdbx_seq_one_letter_code"])
    assert refs == [("PDB", "XXXX", "XXXX", design)]
    seqs = _rows(block, "_struct_ref_seq", ["pdbx_strand_id", "seq_align_beg", "seq_align_end",
                                             "db_align_beg", "db_align_end"])
    assert ("A", "1", str(len(design)), "1", str(len(design))) in seqs
    assert _rows(block, "_struct_ref_seq_dif", ["seq_num"]) == []
    gen = _rows(block, "_entity_src_gen", ["pdbx_gene_src_scientific_name",
                                           "pdbx_gene_src_ncbi_taxonomy_id",
                                           "pdbx_host_org_scientific_name", "pdbx_gene_src_gene"])
    assert gen == [("synthetic construct", "32630", "Escherichia coli", "?")]
    assert _rows(block, "_entity", ["src_method", "details"]) == [("man", spec["design_details"])]
    sw = _rows(block, "_software", ["name", "classification"])
    assert ("RFdiffusion", "model building") in sw and ("REFMAC", "refinement") in sw


def test_a_synthesised_peptide_is_a_synthetic_source(tmp_path):
    peptide = "ACDEFGHIKLMNPQRSTVWY"
    doc, block, st = _prepared(tmp_path, construct=peptide, start=1)
    spec = dict(BASE_SPEC, entities=[{"entity_id": "1", "description": "stapled peptide",
                                      "segments": [{"beg": 1, "end": 20, "kind": "synthetic_peptide"}]}])
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    assert mmcif_edit.split_categories(block) == []
    assert _rows(block, "_entity", ["src_method"]) == [("syn",)]
    assert _rows(block, "_entity_src_gen", ["entity_id"]) == []
    syn = _rows(block, "_pdbx_entity_src_syn", ["entity_id", "pdbx_beg_seq_num", "pdbx_end_seq_num",
                                                "organism_scientific", "ncbi_taxonomy_id"])
    assert syn == [("1", "1", "20", "synthetic construct", "32630")]
    assert _rows(block, "_struct_ref", ["db_name", "pdbx_seq_one_letter_code"]) == [("PDB", peptide)]


def test_a_chimera_of_designed_and_natural_segments(tmp_path):
    design = "SEEALKRLAEEAIKLAKEG"
    construct = design + "GSG" + REFERENCE
    doc, block, st = _prepared(tmp_path, construct=construct, start=1)
    n = len(design)
    spec = dict(BASE_SPEC, entities=[{"entity_id": "1", "description": "binder-target fusion",
                                      "segments": [
                                          {"beg": 1, "end": n, "kind": "designed"},
                                          {"beg": n + 1, "end": n + 3, "kind": "linker"},
                                          {"beg": n + 4, "end": len(construct), "kind": "natural",
                                           "organism": "Thermus thermophilus", "taxid": 274,
                                           "gene": "ilvC", "uniprot": "P00000"}]}])
    # A chimera's classification is the depositor's call: one term, the rest in the text.
    spec["keywords"] = "PROTEIN BINDING"
    spec["keyword_text"] = "de novo binder, fusion"
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    assert mmcif_edit.split_categories(block) == []
    assert gemmi.cif.as_string(block.find_value("_struct_keywords.pdbx_keywords")) == "PROTEIN BINDING"
    gen = _rows(block, "_entity_src_gen", ["pdbx_src_id", "pdbx_seq_type", "pdbx_beg_seq_num",
                                           "pdbx_end_seq_num", "pdbx_gene_src_scientific_name",
                                           "pdbx_gene_src_gene"])
    assert gen == [("1", "Biological sequence", "1", str(n), "synthetic construct", "?"),
                   ("2", "Linker", str(n + 1), str(n + 3), "synthetic construct", "?"),
                   ("3", "Biological sequence", str(n + 4), str(len(construct)),
                    "Thermus thermophilus", "ilvC")]
    # Only the natural segment references UniProt; the designed part is
    # covered by no reference and is not listed as a "difference".
    assert _rows(block, "_struct_ref", ["db_name", "pdbx_db_accession"]) == [("UNP", "P00000")]
    seqs = _rows(block, "_struct_ref_seq", ["pdbx_strand_id", "seq_align_beg", "seq_align_end",
                                             "db_align_beg", "db_align_end"])
    assert ("A", str(n + 4), str(len(construct)), "1", str(len(REFERENCE))) in seqs
    difs = _rows(block, "_struct_ref_seq_dif", ["seq_num", "details"])
    assert all(d[1] == "linker" for d in difs) and len(difs) == 3 * 2  # GSG, both chains


def test_unfetchable_reference_is_recorded_with_a_warning(tmp_path):
    doc, block, st = _prepared(tmp_path)

    def failing(acc):
        raise RuntimeError("UniProt %s: HTTP 500" % acc)

    spec = dict(BASE_SPEC, entities=[{"entity_id": "1", "description": "toy", "segments": [
        {"beg": 9, "end": 33, "kind": "natural", "uniprot": "P00000"}]}])
    warnings = metadata.write_metadata(block, st, spec, fetch=failing)
    assert warnings and "P00000" in warnings[0]
    assert _rows(block, "_struct_ref", ["db_name", "pdbx_db_accession"]) == [("UNP", "P00000")]


def test_statistics_merge_replaces_the_model_tables(tmp_path):
    doc, block, st = _prepared(tmp_path)
    block.set_mmcif_category("_reflns.", {"entry_id": ["XXXX"], "d_resolution_high": ["9.9"]})
    stats = tmp_path / "DataStatistics.cif"
    stats.write_text("data_stats\n_entry.id XXXX\n_cell.length_a 1\n"
                     "_reflns.pdbx_ordinal 1\n_reflns.d_resolution_high 2.10\n"
                     "_reflns.pdbx_CC_half 0.99\n"
                     "loop_\n_reflns_shell.pdbx_ordinal\n_reflns_shell.d_res_high\n"
                     "_reflns_shell.pdbx_CC_half\n1 5.0 0.999\n2 2.1 0.50\n")
    written = statistics.merge_statistics(block, stats)
    assert written == ["_reflns", "_reflns_shell"]
    assert gemmi.cif.as_string(block.find_value("_reflns.d_resolution_high")) == "2.10"
    assert len(block.find_mmcif_category("_reflns_shell.")) == 1   # the outermost shell only
    assert gemmi.cif.as_string(block.find_value("_cell.length_a")) == "50"  # the model's, untouched
    assert mmcif_edit.split_categories(block) == []
    xml = statistics.cifstats_to_xml(stats)
    labels = {s.get("label"): s for s in xml.findall("Stat")}
    assert labels["Half-set correlation CC(1/2)"].get("outer") == "0.50"


def test_reference_alignment_and_fasta_parsing():
    ref = Reference("1", "UNP", "TOY", "P0", REFERENCE, CONSTRUCT, 9, 33, "D8H")
    assert ref.segments == [(9, 22, 1, 14), (26, 33, 15, 22)]
    details = {d[4] for d in ref.diffs}
    assert {DIF_CONFLICT, DIF_ENGINEERED, DIF_INSERTION} <= details
    assert parse_fasta(">sp|Q5SJ03|ILVC_THET8 Ketol-acid\nMKIYY\nEHD\n", "Q5SJ03") == ("ILVC_THET8", "MKIYYEHD")
