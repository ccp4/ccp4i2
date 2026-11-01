"""Findings of the independent review of the first prepare_deposition commits,
each pinned by a test.
"""
import gemmi
import pytest

from ccp4i2.lib.deposition import metadata, mmcif_edit, struct_ref as sr
from ccp4i2.lib.deposition.entities import rebuild_entities_in
from ccp4i2.tests.unit.plugins.test_prepare_deposition_metadata import (
    BASE_SPEC,
    CONSTRUCT,
    _fetch,
    _rows,
    _toy_model,
)

# Thermus thermophilus HB8 ketol-acid reductoisomerase, UniProt Q5SJ03.
Q5SJ03 = (
    "MKIYYEHDADLGFILGKKVAVLGFGSQGHAHALNLKDSGVDVRVGLRKGSRSWEKAEAAG"
    "LRVLPVAEAVREADVVMVLLPDEKQAQVYREEVEPNLKEGGALAFAHGFNVHFGQIKPRK"
    "DLDVWMVAPKGPGHLVRSEYQKGSGVPALVAVHQDASGSAFPTALAYAKAIGAARAGVIA"
    "TTFKDETETDLFGEQAVLCGGLTRLIRAGFETLVEAGYPPEMAYFETVHEVKLIVDLIYE"
    "AGLKGMRYSISNTAEYGDYTRGDLAVPLEETKRRMREILRQIQSGEFAREWMLENQVGSP"
    "VLEANRKRWAAHPIEEVGSRLRAMMPFLKARVMEEVG")
# The Ahir3252 construct: His tag, DARPin H10-2-G3, GGSGGSGG linker, AHIR with
# L33H and a StrepTag-II inserted in place of Q141-K142, RS cloning remnant.
AHIR3252 = (
    "MTGTHHHHHHSDLGKKLLEAARAGQDDEVRILMANGADVNAKDEYGLTPLYLATAHGHLEIVEVLLKNGADVNAVDAIGFTPLHLAAFIG"
    "HLEIAEVLLKHGADVNAQDKFGKTAFDISIGNGNEDLAEILQKLGGSGGSGGLQIYYEHDADLGFILGKKVAVLGFGSQGHAHAHNLKDS"
    "GVDVRVGLRKGSRSWEKAEAAGLRVLPVAEAVREADVVMVLLPDEKQAQVYREEVEPNLKEGGALAFAHGFNVHFGQIKPRKDLDVWMVA"
    "PKGPGHLVRSEYGGGSNWSHPQFEKRPPGGSGVPALVAVHQDASGSAFPTALAYAKAIGAARAGVIATTFKDETETDLFGEQAVLCGGLT"
    "RLIRAGFETLVEAGYPPEMAYFETVHEVKLIVDLIYEAGLKGMRYSISNTAEYGDYTRGDLAVPLEETKRRMREILRQIQSGEFAREWML"
    "ENQVGSPVLEANRKRWAAHPIEEVGSRLRAMMRS")


def test_ahir3252_alignment_segments_and_differences():
    assert len(AHIR3252) == 484 and len(Q5SJ03) == 337
    ref = sr.Reference("1", "UNP", "ILVC_THET8", "Q5SJ03", Q5SJ03, AHIR3252, 143, 484,
                       "L33H, Q141P, K142G")
    # The reference end gaps are free: no bogus deletions past M325, and the
    # tag inserted in place of Q141-K142 is one insertion plus two
    # substitutions, so two mapped stretches.
    assert ref.segments == [(145, 282, 3, 140), (300, 482, 143, 325)]
    details = {}
    for cpos, cres, rpos, rres, det in ref.diffs:
        details[cpos] = (det, rpos)
    assert len(ref.diffs) == 22
    assert details[143] == (sr.DIF_CLONING, None) and details[144] == (sr.DIF_CLONING, None)
    assert details[175] == (sr.DIF_ENGINEERED, 33)
    assert all(details[p] == (sr.DIF_INSERTION, None) for p in range(283, 298))
    assert details[298] == (sr.DIF_ENGINEERED, 141) and details[299] == (sr.DIF_ENGINEERED, 142)
    assert details[483] == (sr.DIF_CLONING, None) and details[484] == (sr.DIF_CLONING, None)
    assert not [d for d in ref.diffs if d[4] == sr.DIF_DELETION]


def test_a_domain_cut_from_the_middle_is_not_a_run_of_deletions():
    piece = "M" + Q5SJ03[100:250]
    ref = sr.Reference("1", "UNP", "X", "X", Q5SJ03, piece, 1, len(piece))
    assert ref.segments == [(2, 151, 101, 250)]
    assert [d[4] for d in ref.diffs] == [sr.DIF_CLONING]


def test_a_his_tag_is_not_a_run_of_conflicts():
    construct = "MGHHHHHH" + Q5SJ03[:60]
    ref = sr.Reference("1", "UNP", "X", "X", Q5SJ03, construct, 1, len(construct))
    assert ref.segments == [(9, 68, 1, 60)]
    assert {d[4] for d in ref.diffs} == {sr.DIF_CLONING} and len(ref.diffs) == 8


def test_author_numbering_is_written_not_construct_numbering(tmp_path):
    # chains modelled from construct residue 3; author numbering shifted by +10
    path = _toy_model(tmp_path / "model.cif", CONSTRUCT)
    doc = gemmi.cif.read(str(path))
    st0 = gemmi.make_structure_from_block(doc[0])
    for chain in st0[0]:
        for res in chain:
            res.seqid = gemmi.SeqId(res.seqid.num + 10, " ")
    doc = st0.make_mmcif_document()
    doc[0].set_mmcif_category("_diffrn_radiation_wavelength.", {"id": ["1"], "wavelength": ["0.9795"]})
    st = rebuild_entities_in(doc, [("construct", CONSTRUCT, "PROTEIN")])
    block = doc[0]
    auth = sr.auth_numbering(st)
    assert auth["A"][3] == 13
    spec = dict(BASE_SPEC, entities=[{"entity_id": "1", "description": "toy", "segments": [
        {"beg": 9, "end": 33, "kind": "engineered", "uniprot": "P00000", "engineered": "D8H"}]}])
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    seqs = _rows(block, "_struct_ref_seq", ["pdbx_strand_id", "seq_align_beg",
                                             "pdbx_auth_seq_align_beg", "pdbx_auth_seq_align_end"])
    row = next(r for r in seqs if r[0] == "A")
    assert int(row[2]) == int(row[1]) + 10
    difs = _rows(block, "_struct_ref_seq_dif", ["pdbx_pdb_strand_id", "seq_num", "pdbx_auth_seq_num"])
    assert all(int(r[2]) == int(r[1]) + 10 for r in difs if r[0] == "A")


def test_wavelength_from_the_form_fills_a_model_without_one(tmp_path):
    # A PDB-format model carries no wavelength; the form's is written and
    # only its absence after that is a problem.
    path = _toy_model(tmp_path / "model.cif", CONSTRUCT)
    doc = gemmi.cif.read(str(path))
    doc[0].find_mmcif_category("_diffrn_radiation_wavelength.").erase()
    st = rebuild_entities_in(doc, [("construct", CONSTRUCT, "PROTEIN")])
    spec = dict(BASE_SPEC, experiment=dict(BASE_SPEC["experiment"], wavelength=0.9795))
    warnings = metadata.write_metadata(doc[0], st, spec, fetch=_fetch)
    assert mmcif_edit.value_of(doc[0], "_diffrn_radiation_wavelength.wavelength") == "0.9795"
    assert not [w for w in warnings if "wavelength" in w]
    doc = gemmi.cif.read(str(path))
    doc[0].find_mmcif_category("_diffrn_radiation_wavelength.").erase()
    st = rebuild_entities_in(doc, [("construct", CONSTRUCT, "PROTEIN")])
    warnings = metadata.write_metadata(doc[0], st, BASE_SPEC, fetch=_fetch)
    assert any("no wavelength" in w for w in warnings)


def test_mad_wavelengths_are_preserved(tmp_path):
    path = _toy_model(tmp_path / "model.cif", CONSTRUCT)
    doc = gemmi.cif.read(str(path))
    doc[0].set_mmcif_category("_diffrn_radiation_wavelength.", {"id": ["1", "2", "3"],
                                                                "wavelength": ["0.9795", "0.9793", "0.9500"]})
    st = rebuild_entities_in(doc, [("construct", CONSTRUCT, "PROTEIN")])
    spec = dict(BASE_SPEC, experiment=dict(BASE_SPEC["experiment"], wavelength=0.9795))
    warnings = metadata.write_metadata(doc[0], st, spec, fetch=_fetch)
    assert len(doc[0].find_mmcif_category("_diffrn_radiation_wavelength.")) == 3
    assert any("3 wavelengths" in w for w in warnings)
    assert mmcif_edit.split_categories(doc[0]) == []


def test_set_pairs_in_place_on_a_one_row_loop_keeps_the_other_items():
    block = gemmi.cif.read_string(
        "data_x\nloop_\n_refine.entry_id\n_refine.ls_d_res_high\n_refine.ls_R_factor_R_work\n"
        "XXXX 2.1 0.19\n_cell.length_a 50\n")[0]
    mmcif_edit.set_pairs_in_place(block, "_refine", [("pdbx_starting_model", "2JAB")])
    tags, rows = mmcif_edit.category_rows(block, "_refine")
    assert tags == ["entry_id", "ls_d_res_high", "ls_R_factor_R_work", "pdbx_starting_model"]
    assert rows == [["XXXX", "2.1", "0.19", "2JAB"]]
    assert mmcif_edit.split_categories(block) == []


def test_a_chain_unlike_every_au_sequence_keeps_its_own_entity(tmp_path):
    path = _toy_model(tmp_path / "model.cif", CONSTRUCT, chains=("A",))
    doc = gemmi.cif.read(str(path))
    st0 = gemmi.make_structure_from_block(doc[0])
    peptide = gemmi.Chain("P")
    for num, one in enumerate("WWPPWWPPWW", start=1):
        res = gemmi.Residue()
        res.name = gemmi.find_tabulated_residue(gemmi.expand_one_letter(one, gemmi.ResidueKind.AA)).name
        res.seqid = gemmi.SeqId(num, " ")
        atom = gemmi.Atom()
        atom.name, atom.element, atom.pos = "CA", gemmi.Element("C"), gemmi.Position(num, 50, 0)
        atom.occ, atom.b_iso = 1.0, 30.0
        res.add_atom(atom)
        peptide.add_residue(res)
    st0[0].add_chain(peptide)
    st0.setup_entities()
    doc = st0.make_mmcif_document()
    warnings = []
    st = rebuild_entities_in(doc, [("construct", CONSTRUCT, "PROTEIN")], warnings=warnings)
    polymers = [e for e in st.entities if e.entity_type == gemmi.EntityType.Polymer]
    assert len(polymers) == 2
    assert warnings and "chain P" in warnings[0]
    strands = {row.str(0): row.str(1) for row in doc[0].find("_entity_poly.", ["entity_id", "pdbx_strand_id"])}
    assert strands["1"] == "A" and "P" in strands.values()


def test_crystallisation_method_is_free_text_because_the_values_hold_commas():
    from ccp4i2.core.tasks import locate_def_xml
    from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser

    container = DefXmlParser().parse_def_xml(locate_def_xml("prepare_deposition"))
    ctrl = container.controlParameters
    assert not ctrl.CRYSTAL_METHOD.get_qualifier("enumerators")
    assert bool(ctrl.KEEP_HYDROGENS) is True   # the default: deposit the refined model as validated
