"""What the DepUI showed on the first deposition made from a metadata-bearing
file (2026-11-01), fixed and pinned:

- a ``_software.name`` outside OneDep's list is dropped silently;
- the starting model table is ``_pdbx_initial_refinement_model``;
- OneDep reads the ``_obs`` statistics items and the outermost shell;
- a ligand of interest is declared in ``_pdbx_entity_instance_feature``;
- the assembly's experimental support is ``_pdbx_struct_assembly_auth_evidence``;
- ``_struct_keywords.pdbx_keywords`` is one term of the classification list.
"""
import gemmi
import pytest

from ccp4i2.lib.deposition import metadata, mmcif_edit, statistics, vocab
from ccp4i2.lib.deposition.entities import rebuild_entities_in
from ccp4i2.tests.unit.plugins.test_prepare_deposition_metadata import (
    BASE_SPEC,
    CONSTRUCT,
    _fetch,
    _rows,
    _toy_model,
)


def _prepared_with_ligand(tmp_path):
    path = _toy_model(tmp_path / "model.cif", CONSTRUCT)
    doc = gemmi.cif.read(str(path))
    st = gemmi.make_structure_from_block(doc[0])
    # two chloride ions, one per chain, as HETATM residues
    for chain in st[0]:
        res = gemmi.Residue()
        res.name, res.seqid, res.het_flag = "CL", gemmi.SeqId(201, " "), "H"
        atom = gemmi.Atom()
        atom.name, atom.element, atom.pos = "CL", gemmi.Element("Cl"), gemmi.Position(1, 2, 3)
        atom.occ, atom.b_iso = 1.0, 30.0
        res.add_atom(atom)
        chain.add_residue(res)
    st.setup_entities()
    doc = st.make_mmcif_document()
    doc[0].set_mmcif_category("_diffrn_radiation_wavelength.", {"id": ["1"], "wavelength": ["0.9795"]})
    doc[0].set_mmcif_category("_refine.", {"entry_id": ["XXXX"], "ls_d_res_high": ["2.1"]})
    st = rebuild_entities_in(doc, [("construct", CONSTRUCT, "PROTEIN")])
    return doc[0], st


def test_software_names_follow_onedep_spelling_and_the_rest_is_reported(tmp_path):
    block, st = _prepared_with_ligand(tmp_path)
    spec = dict(BASE_SPEC, software=[
        {"name": "refmac", "version": "5.8.0431", "classification": "refinement"},
        {"name": "CCP4i2", "version": "3.1.0a109", "classification": "data processing"},
        {"name": "aimless", "version": "0.8.2", "classification": "data scaling"},
        {"name": "my_script.py", "classification": "other"},
        {"name": "ProteinMPNN", "classification": "model building"}])
    warnings = metadata.write_metadata(block, st, spec, fetch=_fetch)
    names = [r[0] for r in _rows(block, "_software", ["name"])]
    # Unverified names are written (the archive records what was used), with a warning.
    assert names == ["REFMAC", "Aimless", "my_script.py", "ProteinMPNN"]
    assert any("CCP4i2" in w and "not written" in w for w in warnings)
    assert any("my_script.py" in w and "Communication tab" in w for w in warnings)
    assert any("ProteinMPNN" in w for w in warnings)
    assert vocab.canonical_software_name("Refmacat") == "REFMAC"
    assert vocab.canonical_software_name("ccp4i2") is None


def test_starting_models_table_and_summary(tmp_path):
    block, st = _prepared_with_ligand(tmp_path)
    spec = dict(BASE_SPEC, phasing={"method": "MR", "starting_models": [
        {"type": "experimental model", "source": "PDB", "accession": "2JAB"},
        {"type": "experimental model", "source": "Other", "details": "in-house AHIR dimer"},
        {"type": "in silico model", "source": "AlphaFold", "accession": "AF-Q5SJ03-F1", "entity_ids": "1"}]})
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    rows = _rows(block, "_pdbx_initial_refinement_model",
                 ["id", "entity_id_list", "type", "source_name", "accession_code", "details"])
    assert rows == [("1", "1", "experimental model", "PDB", "2JAB", "?"),
                    ("2", "1", "experimental model", "Other", "?", "in-house AHIR dimer"),
                    ("3", "1", "in silico model", "AlphaFold", "AF-Q5SJ03-F1", "?")]
    assert gemmi.cif.as_string(block.find_value("_refine.pdbx_starting_model")) == \
        "2JAB; in-house AHIR dimer; AF-Q5SJ03-F1"
    assert mmcif_edit.split_categories(block) == []


def test_obs_items_and_outer_shell(tmp_path):
    stats = tmp_path / "DataStatistics.cif"
    stats.write_text(
        "data_stats\n_reflns.pdbx_ordinal 1\n_reflns.d_resolution_high 2.10\n"
        "_reflns.number_obs 115300\n_reflns.pdbx_Rmerge_I_all 0.409\n"
        "_reflns.pdbx_netI_over_sigmaI 4.5\n_reflns.percent_possible_obs 91.4\n"
        "loop_\n_reflns_shell.pdbx_ordinal\n_reflns_shell.d_res_low\n_reflns_shell.d_res_high\n"
        "_reflns_shell.number_unique_all\n_reflns_shell.Rmerge_I_all\n"
        "_reflns_shell.pdbx_netI_over_sigmaI_all\n_reflns_shell.percent_possible_all\n"
        "_reflns_shell.pdbx_CC_half\n"
        "1 81.56 14.07 768 0.062 16.5 85.8 0.997\n"
        "2 14.07 9.93 1313 0.067 18.6 86.8 0.998\n"
        "3 2.19 2.10 7609 1.90 0.9 88.1 0.41\n")
    doc = gemmi.cif.read_string("data_x\n_entry.id XXXX\n")
    written = statistics.merge_statistics(doc[0], stats)
    assert written == ["_reflns", "_reflns_shell"]
    shell = _rows(doc[0], "_reflns_shell", ["pdbx_ordinal", "d_res_high", "number_unique_obs",
                                            "Rmerge_I_obs", "meanI_over_sigI_obs",
                                            "percent_possible_obs", "number_unique_all"])
    assert shell == [("1", "2.10", "7609", "1.90", "0.9", "88.1", "7609")]
    assert gemmi.cif.as_string(doc[0].find_value("_reflns.pdbx_Rmerge_I_obs")) == "0.409"
    # all shells on request, every row with its _obs twins
    doc = gemmi.cif.read_string("data_x\n_entry.id XXXX\n")
    statistics.merge_statistics(doc[0], stats, all_shells=True)
    assert len(doc[0].find_mmcif_category("_reflns_shell.")) == 3
    assert mmcif_edit.split_categories(doc[0]) == []


def test_ligand_of_interest_is_declared_per_instance(tmp_path):
    block, st = _prepared_with_ligand(tmp_path)
    spec = dict(BASE_SPEC, ligand_of_interest="cl")
    warnings = metadata.write_metadata(block, st, spec, fetch=_fetch)
    rows = _rows(block, "_pdbx_entity_instance_feature",
                 ["ordinal", "feature_type", "auth_asym_id", "auth_seq_num", "comp_id"])
    assert rows == [("1", "SUBJECT OF INVESTIGATION", "A", "201", "CL"),
                    ("2", "SUBJECT OF INVESTIGATION", "B", "201", "CL")]
    assert not [w for w in warnings if "ligand" in w]
    # none declared: nothing written (the report tells the user to tick the box)
    block, st = _prepared_with_ligand(tmp_path)
    metadata.write_metadata(block, st, BASE_SPEC, fetch=_fetch)
    assert _rows(block, "_pdbx_entity_instance_feature", ["ordinal"]) == []
    # a code not in the model is a warning
    block, st = _prepared_with_ligand(tmp_path)
    warnings = metadata.write_metadata(block, st, dict(BASE_SPEC, ligand_of_interest="ATP"), fetch=_fetch)
    assert any("ATP" in w for w in warnings)


def test_assembly_evidence_table(tmp_path):
    block, st = _prepared_with_ligand(tmp_path)
    metadata.write_metadata(block, st, BASE_SPEC, fetch=_fetch)
    assert _rows(block, "_pdbx_struct_assembly_auth_evidence",
                 ["id", "assembly_id", "experimental_support", "details"]) == [("1", "1", "none", "?")]
    assert _rows(block, "_pdbx_struct_assembly", ["id", "details"]) == [("1", "author_defined_assembly")]
    block, st = _prepared_with_ligand(tmp_path)
    spec = dict(BASE_SPEC, assembly={"mode": "author", "oligomeric_count": 2, "operators": "",
                                     "evidence": "gel filtration", "details": "dimer on S200"})
    metadata.write_metadata(block, st, spec, fetch=_fetch)
    assert _rows(block, "_pdbx_struct_assembly_auth_evidence",
                 ["experimental_support", "details"]) == [("gel filtration", "dimer on S200")]
    block, st = _prepared_with_ligand(tmp_path)
    warnings = metadata.write_metadata(block, st, dict(BASE_SPEC, assembly={"evidence": "PISA"}), fetch=_fetch)
    assert any("PISA" in w for w in warnings)
    assert _rows(block, "_pdbx_struct_assembly_auth_evidence", ["experimental_support"]) == [("none",)]


def test_classification_is_one_enumerated_term(tmp_path):
    block, st = _prepared_with_ligand(tmp_path)
    with pytest.raises(ValueError, match="one term"):
        metadata.write_metadata(block, st, dict(BASE_SPEC, keywords="OXIDOREDUCTASE/DE NOVO PROTEIN"),
                                fetch=_fetch)
    block, st = _prepared_with_ligand(tmp_path)
    with pytest.raises(ValueError, match="classification"):
        metadata.write_metadata(block, st, dict(BASE_SPEC, keywords=None), fetch=_fetch)
    block, st = _prepared_with_ligand(tmp_path)
    metadata.write_metadata(block, st, dict(BASE_SPEC, keywords="oxidoreductase",
                                            keyword_text="DE NOVO PROTEIN, DARPin fusion"), fetch=_fetch)
    assert gemmi.cif.as_string(block.find_value("_struct_keywords.pdbx_keywords")) == "OXIDOREDUCTASE"
    assert "DARPin" in gemmi.cif.as_string(block.find_value("_struct_keywords.text"))
    # wholly designed: DE NOVO PROTEIN without being asked; a chimera is not wholly designed
    assert metadata.default_keywords([{"segments": [{"kind": "tag"}, {"kind": "designed"}]}]) == "DE NOVO PROTEIN"
    assert metadata.default_keywords([{"segments": [{"kind": "designed"}, {"kind": "natural"}]}]) is None


def test_the_def_xml_lists_are_the_vocabulary():
    from ccp4i2.core.tasks import locate_def_xml
    from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser

    container = DefXmlParser().parse_def_xml(locate_def_xml("prepare_deposition"))
    ctrl = container.controlParameters
    # The def.xml holds no lists for these: the plugin gives them their menus
    # from the one source (see test_prepare_deposition_vocabulary).
    assert not ctrl.KEYWORDS.get_qualifier("enumerators")
    assert not ctrl.ASSEMBLY_EVIDENCE.get_qualifier("enumerators")
    sub = ctrl.STARTING_MODELS.get_qualifier("subItem")
    assert sub["class"].__name__ == "CStartingModel"
