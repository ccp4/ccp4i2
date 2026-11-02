"""What the sanitising stage does to the mmCIF each refinement program writes.

One fixture per program, written as the program writes it (as far as a
fixture can): Refmac/refmacat with its private categories, dataset-named
wavelength and appended dictionary blocks; servalcat with per-shell R only;
BUSTER without ``pdbx_refine_id`` (defensive: no BUSTER mmCIF was reachable).
Every fix comes back as "found X, did Y" and the tests read those.
"""
import json

import gemmi
import numpy
import pytest

from ccp4i2.lib.deposition import sanitise

ATOMS = """
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 1 1 ? 0.0 0.0 0.0 1.00 20.0 1 ALA A N 1
ATOM 2 C CA . ALA A 1 1 ? 1.5 0.0 0.0 1.00 20.0 1 ALA A CA 1
ATOM 3 C C . ALA A 1 1 ? 2.0 1.4 0.0 1.00 20.0 1 ALA A C 1
ATOM 4 O O . ALA A 1 1 ? 1.2 2.3 0.0 1.00 20.0 1 ALA A O 1
ATOM 5 H H . ALA A 1 1 ? -0.5 0.5 0.5 1.00 20.0 1 ALA A H 1
ATOM 6 N N . GLY A 1 2 ? 3.3 1.5 0.0 {occ} 20.0 2 GLY A N 1
ATOM 7 C CA . GLY A 1 2 ? 4.0 2.8 0.0 1.00 20.0 2 GLY A CA 1
ATOM 8 C C . GLY A 1 2 ? 5.5 2.7 0.0 1.00 20.0 2 GLY A C 1
ATOM 9 O O . GLY A 1 2 ? 6.1 1.6 0.0 1.00 20.0 2 GLY A O 1
"""

HEADER = """data_old
_entry.id 1ABC
_cell.length_a 50.0
_cell.length_b 60.0
_cell.length_c 70.0
_cell.angle_alpha 90.0
_cell.angle_beta 90.0
_cell.angle_gamma 90.0
_symmetry.space_group_name_H-M 'P 1'
_exptl.entry_id 1ABC
_exptl.method 'X-RAY DIFFRACTION'
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_strand_id
1 polypeptide(L) A
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
_pdbx_database_status.entry_id 1ABC
_pdbx_database_status.recvd_initial_deposition_date 2015-12-22
loop_
_audit_author.name
_audit_author.pdbx_ordinal
'Someone, E.' 1
'Else, A.N.' 2
_citation.id primary
_citation.title 'An earlier structure the search model came from'
"""

REFMAC = HEADER + """
loop_
_software.pdbx_ordinal
_software.name
_software.version
_software.classification
1 refmac 5.8.0431 refinement
_refine.entry_id 1ABC
_refine.ls_d_res_high 2.1
_refine.ls_R_factor_R_work 0.19
_refine.ls_R_factor_R_free 0.23
_ccp4_refine_ls.jelly_body 1
_ccp4_refine_tls.details 'added to B'
_ccp4_struct_ncs.method local
loop_
_diffrn_radiation_wavelength.id
_diffrn_radiation_wavelength.wavelength
FROMDIALS 0.9795
_diffrn_radiation.diffrn_id 1
_diffrn_radiation.wavelength_id FROMDIALS
""" + ATOMS.format(occ="1.00") + """
data_comp_list
_chem_comp.id ALA
data_ALA
_chem_comp_atom.comp_id ALA
"""

SERVALCAT = HEADER + """
loop_
_software.pdbx_ordinal
_software.name
_software.version
_software.classification
1 servalcat 0.4.123 refinement
_refine.entry_id 1ABC
_refine.pdbx_refine_id 'X-RAY DIFFRACTION'
_refine.ls_d_res_high 2.56
loop_
_refine_ls_shell.pdbx_refine_id
_refine_ls_shell.d_res_high
_refine_ls_shell.number_reflns_R_work
_refine_ls_shell.number_reflns_R_free
_refine_ls_shell.R_factor_R_work
_refine_ls_shell.R_factor_R_free
'X-RAY DIFFRACTION' 5.0 1000 50 0.15 0.18
'X-RAY DIFFRACTION' 2.56 2000 100 0.25 0.30
_diffrn_radiation_wavelength.id 1
_diffrn_radiation_wavelength.wavelength 0.9795
""" + ATOMS.format(occ="1.00")

BUSTER = HEADER + """
loop_
_software.pdbx_ordinal
_software.name
_software.version
_software.classification
1 BUSTER 2.11.6 refinement
_refine.entry_id 1ABC
_refine.ls_d_res_high 2.1
_refine.ls_R_factor_R_work 0.20
_refine.ls_R_factor_R_free 0.24
_refine.details 'BUSTER-TNT with LSSR restraints'
_diffrn_radiation_wavelength.id 1
_diffrn_radiation_wavelength.wavelength 0.9795
""" + ATOMS.format(occ="1.00")


def _found(fixes):
    return "\n".join("%s | %s" % (f.found, f.did) for f in fixes)


def test_refmac_private_categories_blocks_and_wavelength_id():
    doc = gemmi.cif.read_string(REFMAC)
    assert len(doc) == 3
    program, fixes, problems, st = sanitise.sanitise(doc)
    text = _found(fixes)
    assert program == "refmac" and problems == []
    assert len(doc) == 1
    block = doc[0]
    for cat in sanitise.REFMAC_PRIVATE:
        assert not block.find_mmcif_category(cat + "."), cat
    assert "3 data blocks" in text and "kept only the model block" in text
    assert "_ccp4_refine_ls" in text
    assert gemmi.cif.as_string(block.find_value("_diffrn_radiation_wavelength.id")) == "1"
    assert gemmi.cif.as_string(block.find_value("_diffrn_radiation.wavelength_id")) == "1"
    assert gemmi.cif.as_string(block.find_value("_diffrn_radiation_wavelength.wavelength")) == "0.9795"
    assert "FROMDIALS" in text


def test_inherited_header_is_stripped_and_named():
    doc = gemmi.cif.read_string(REFMAC)
    _program, fixes, _problems, _st = sanitise.sanitise(doc)
    block = doc[0]
    text = _found(fixes)
    for cat in ("_pdbx_database_status", "_audit_author", "_citation"):
        assert not block.find_mmcif_category(cat + "."), cat
        assert cat in text
    assert "Someone, E." in text          # the stale author list is named, not lost silently
    assert "earlier structure" in text    # and so is the stale citation
    assert gemmi.cif.as_string(block.find_value("_entry.id")) == "XXXX"
    assert "_entry.id 1ABC" in text


def test_hydrogens_are_kept_unless_asked_to_strip():
    # The deposited model is the refined and validated one: hydrogens stay.
    doc = gemmi.cif.read_string(REFMAC)
    _p, fixes, _problems, st = sanitise.sanitise(doc)
    assert "kept" in _found(fixes) and st[0].count_atom_sites() == 9
    doc = gemmi.cif.read_string(REFMAC)
    _p, fixes, _problems, st = sanitise.sanitise(doc, keep_hydrogens=False)
    assert "computed with them present" in _found(fixes)
    assert st[0].count_atom_sites() == 8


def test_bad_occupancy_is_a_problem_not_a_fix():
    doc = gemmi.cif.read_string(HEADER + "_refine.entry_id 1ABC\n"
                                "_diffrn_radiation_wavelength.id 1\n"
                                "_diffrn_radiation_wavelength.wavelength 0.98\n"
                                + ATOMS.format(occ="1.30"))
    _p, _fixes, problems, _st = sanitise.sanitise(doc)
    assert any("occupancy" in p for p in problems)


def test_missing_wavelength_is_not_a_problem_at_this_stage():
    # A PDB-format model never carries one; the task writes it from the form
    # or the reflection data and checks afterwards.
    doc = gemmi.cif.read_string(HEADER + "_refine.entry_id 1ABC\n" + ATOMS.format(occ="1.00"))
    _p, _fixes, problems, _st = sanitise.sanitise(doc)
    assert not any("wavelength" in p for p in problems)


def _servalcat_job(tmp_path, summary):
    job = tmp_path / "job_6"
    sub = job / "job_6_1"
    sub.mkdir(parents=True)
    (sub / "refined_stats.json").write_text(json.dumps([
        {"data": {"summary": {"Rwork": 0.99, "Rfree": 0.99}}},   # an earlier cycle
        {"data": {"summary": summary}}]))
    return job


def test_servalcat_overall_r_is_harvested_when_absent(tmp_path):
    job = _servalcat_job(tmp_path, {"Rwork": 0.1902, "Rfree": 0.2297, "R1work": 0.17})
    doc = gemmi.cif.read_string(SERVALCAT)
    program, fixes, problems, _st = sanitise.sanitise(doc, refinement_job_dir=str(job))
    block = doc[0]
    assert program == "servalcat" and problems == []
    assert gemmi.cif.as_string(block.find_value("_refine.ls_R_factor_R_work")) == "0.1902"
    assert gemmi.cif.as_string(block.find_value("_refine.ls_R_factor_R_free")) == "0.2297"
    assert gemmi.cif.as_string(block.find_value("_refine.ls_number_reflns_R_work")) == "3000"
    assert gemmi.cif.as_string(block.find_value("_refine.ls_number_reflns_R_free")) == "150"
    assert gemmi.cif.as_string(block.find_value("_refine.ls_number_reflns_obs")) == "3150"
    text = _found(fixes)
    assert "Rwork/Rfree of refined_stats.json" in text   # the conventional pair was preferred
    assert "R1" not in text.split("Rwork/Rfree")[0]
    # _refine stays one contiguous category after the insertion.
    from ccp4i2.lib.deposition.mmcif_edit import split_categories
    assert split_categories(block) == []


def test_servalcat_r1_goes_to_the_intensity_r_item_not_r_work(tmp_path):
    job = _servalcat_job(tmp_path, {"R1work": 0.17, "R1free": 0.21})
    doc = gemmi.cif.read_string(SERVALCAT)
    _p, fixes, problems, _st = sanitise.sanitise(doc, refinement_job_dir=str(job))
    assert problems == []
    assert doc[0].find_value("_refine.ls_R_factor_R_work") is None
    assert gemmi.cif.as_string(doc[0].find_value("_refine.ls_R_I_factor_obs")) == "0.1700"
    assert "R1" in _found(fixes) and "NotAvailable" in _found(fixes)


def test_servalcat_stats_come_from_the_job_not_the_newest_file_anywhere(tmp_path):
    job = _servalcat_job(tmp_path, {"Rwork": 0.2000, "Rfree": 0.2400})
    deeper = job / "job_6_1" / "scratch"
    deeper.mkdir()
    (deeper / "refined_stats.json").write_text(json.dumps([{"data": {"summary": {"Rwork": 0.01, "Rfree": 0.02}}}]))
    doc = gemmi.cif.read_string(SERVALCAT)
    sanitise.sanitise(doc, refinement_job_dir=str(job))
    assert gemmi.cif.as_string(doc[0].find_value("_refine.ls_R_factor_R_work")) == "0.2000"


def test_servalcat_without_any_r_is_refused(tmp_path):
    doc = gemmi.cif.read_string(SERVALCAT)
    _p, _fixes, problems, _st = sanitise.sanitise(doc, refinement_job_dir=str(tmp_path))
    assert problems and "R-work" in problems[0]


def test_servalcat_that_wrote_its_own_r_is_left_alone(tmp_path):
    text = SERVALCAT.replace("_refine.ls_d_res_high 2.56",
                             "_refine.ls_d_res_high 2.56\n_refine.ls_R_factor_R_work 0.2101\n"
                             "_refine.ls_R_factor_R_free 0.2503")
    job = _servalcat_job(tmp_path, {"Rwork": 0.1, "Rfree": 0.1})   # must NOT be used
    doc = gemmi.cif.read_string(text)
    _p, fixes, problems, _st = sanitise.sanitise(doc, refinement_job_dir=str(job))
    assert problems == []
    assert gemmi.cif.as_string(doc[0].find_value("_refine.ls_R_factor_R_work")) == "0.2101"
    assert "left as written" in _found(fixes)


def test_buster_gets_a_refine_id():
    doc = gemmi.cif.read_string(BUSTER)
    program, fixes, problems, _st = sanitise.sanitise(doc)
    assert program == "buster" and problems == []
    assert gemmi.cif.as_string(doc[0].find_value("_refine.pdbx_refine_id")) == "X-RAY DIFFRACTION"
    assert "pdbx_refine_id" in _found(fixes)


def test_cell_check_against_reflections(tmp_path):
    doc = gemmi.cif.read_string(BUSTER)
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 1")
    mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))
    mtz.add_dataset("x")
    mtz.add_column("F", "F")
    mtz.set_data(numpy.array([[0, 0, 1, 10.0]], dtype=numpy.float32))
    path = tmp_path / "x.mtz"
    mtz.write_to_file(str(path))
    assert sanitise.cell_matches_reflections(doc[0], path) is None
    mtz.set_cell_for_all(gemmi.UnitCell(51, 60, 70, 90, 90, 90))
    mtz.write_to_file(str(path))
    assert "cell a" in sanitise.cell_matches_reflections(doc[0], path)


def test_reflections_cif_is_written_from_the_mini_mtz_labels(tmp_path):
    from ccp4i2.lib.deposition.reflections import write_reflections_cif

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 1")
    mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))
    mtz.add_dataset("FROMDIALS")
    for label, kind in (("F_SIGF_F", "F"), ("F_SIGF_SIGF", "Q"), ("FREERFLAG_FREER", "I")):
        mtz.add_column(label, kind)
    mtz.set_data(numpy.array([[0, 0, 1, 10.0, 1.0, 0], [0, 0, 2, 12.0, 1.1, 1]], dtype=numpy.float32))
    path = tmp_path / "merged.mtz"
    mtz.write_to_file(str(path))
    out = tmp_path / "refl.cif"
    renamed = write_reflections_cif(path, out, log=lambda *a: None)
    assert any("FREERFLAG_FREER -> FREE" in r for r in renamed)
    block = gemmi.cif.read(str(out))[0]
    tags = [t for t in block.find_loop_item("_refln.index_h").loop.tags] if block.find_loop_item("_refln.index_h") else []
    assert "_refln.F_meas_au" in tags and "_refln.status" in tags   # gemmi writes the free set as status o/f
    assert block.find_loop_item("_refln.index_h").loop.length() == 2
