"""Staging produces the contract's input tree, exactly (design note §3.2,
assertions in §15.2): ``datasets/xtal-NNNN/`` with clean integer names,
``Projects.csv`` present and parseable, a relabelled FreeR column where the
source needed one, a re-spelled dictionary where that did; and a manifest
keyed on uuid. Diffable against the real 201-dataset tree when the volume
holding it is mounted. CCP4-free: gemmi and file operations only.
"""
import filecmp
import json
import os
from pathlib import Path

import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.wrappers.pandda_campaign.script.pandda_staging import (
    DATASETS_DIR, DICT_NAME, MANIFEST_JSON, MODEL_NAME, PROJECTS_CSV,
    REFLECTIONS_NAME, DatasetSpec, link_or_copy, read_manifest,
    read_projects_csv, stage_datasets, xtal_name)
from .conftest import LABELS, pdbx_spelling, source_files

REAL_TREE = Path("/Volumes/LocalStore/pandda/BAZ2B")


def _tree(root):
    root = Path(root)
    return sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file())


def test_staged_tree_matches_the_contract(tmp_path, mini_specs):
    root = tmp_path / "staging"
    manifest = stage_datasets(mini_specs, root)
    expected = [PROJECTS_CSV, MANIFEST_JSON]
    for i in range(3):
        for name in (DICT_NAME, REFLECTIONS_NAME, MODEL_NAME):
            expected.append(f"{DATASETS_DIR}/{xtal_name(i)}/{name}")
    assert _tree(root) == sorted(expected)
    # clean integer names: PanDDA takes the crystal number from the last digits
    for entry in manifest["datasets"]:
        digits = entry["xtal"].rsplit("-", 1)[1]
        assert digits.isdigit() and int(digits) < 100000, entry["xtal"]


def test_xtal_numbering_follows_submission_order(tmp_path, mini_specs):
    reordered = list(reversed(mini_specs))
    manifest = stage_datasets(reordered, tmp_path / "staging")
    assert [e["label"] for e in manifest["datasets"]] == [s.label for s in reordered]
    assert [e["xtal"] for e in manifest["datasets"]] == [xtal_name(i) for i in range(3)]


def test_projects_csv_is_what_the_export_writes(tmp_path, mini_specs):
    root = tmp_path / "staging"
    stage_datasets(mini_specs, root)
    text = (root / PROJECTS_CSV).read_text()
    assert text.splitlines()[0] == "Dataset, Project"
    assert read_projects_csv(root / PROJECTS_CSV) == [(xtal_name(i), LABELS[i]) for i in range(3)]


def test_staged_bytes_equal_the_sources(tmp_path, mini_specs):
    root = tmp_path / "staging"
    manifest = stage_datasets(mini_specs, root)
    for spec, entry in zip(mini_specs, manifest["datasets"]):
        dataset_dir = root / DATASETS_DIR / entry["xtal"]
        assert filecmp.cmp(spec.xyzin, dataset_dir / MODEL_NAME, shallow=False)
        # these sources already carry FreeR_flag and type: nothing rewritten
        assert filecmp.cmp(spec.hklin, dataset_dir / REFLECTIONS_NAME, shallow=False)
        assert filecmp.cmp(spec.dictionary, dataset_dir / DICT_NAME, shallow=False)
        assert {f["how"] for f in entry["files"].values()} <= {"link", "copy"}
        assert not any(f["how"] == "rewritten" for f in entry["files"].values())


def test_manifest_is_keyed_on_uuid_and_records_digests(tmp_path, mini_specs):
    root = tmp_path / "staging"
    written = stage_datasets(mini_specs, root, provenance={"contract": "test"})
    manifest = read_manifest(root / MANIFEST_JSON)
    assert manifest == written
    assert manifest["provenance"] == {"contract": "test"}
    for spec, entry in zip(mini_specs, manifest["datasets"]):
        assert entry["project_uuid"] == spec.project_uuid
        for name, info in entry["files"].items():
            staged = root / DATASETS_DIR / entry["xtal"] / name
            assert len(info["sha256"]) == 64
            assert Path(info["source"]).is_file()
            import hashlib
            assert info["sha256"] == hashlib.sha256(staged.read_bytes()).hexdigest()


def test_freer_is_relabelled_where_the_source_needs_it(tmp_path, mini_specs):
    spec = mini_specs[0]
    mtz = gemmi.read_mtz_file(str(spec.hklin))
    col = next(c for c in mtz.columns if c.label == "FreeR_flag")
    col.label = "FREER"                      # CCP4i2's own convention
    freer_src = tmp_path / "freer.mtz"
    mtz.write_to_file(str(freer_src))
    specs = [DatasetSpec(label="needs-relabel", xyzin=spec.xyzin, hklin=freer_src)]
    root = tmp_path / "staging"
    manifest = stage_datasets(specs, root)
    staged = root / DATASETS_DIR / xtal_name(0) / REFLECTIONS_NAME
    labels = [c.label for c in gemmi.read_mtz_file(str(staged)).columns]
    assert "FreeR_flag" in labels and "FREER" not in labels
    entry = manifest["datasets"][0]["files"][REFLECTIONS_NAME]
    assert entry["how"] == "rewritten"
    assert not list((root / DATASETS_DIR / xtal_name(0)).glob("*_pandda.mtz")), "no scratch left behind"
    assert DICT_NAME not in manifest["datasets"][0]["files"], "optional dictionary absent"


def test_freer_is_found_by_column_type_under_ccp4i2s_own_naming(tmp_path, mini_specs):
    """dimple's COMPLETE_MTZ names the free set FREERFLAG_FREER (the parameter
    that supplied it, joined mini-MTZ style): no label list anticipates it,
    and PanDDA accepts only three exact labels. Found by column type."""
    spec = mini_specs[0]
    mtz = gemmi.read_mtz_file(str(spec.hklin))
    for col in mtz.columns:
        if col.label == "FreeR_flag":
            col.label = "FREERFLAG_FREER"
        elif col.label == "F":
            col.label = "F_SIGF_F"
        elif col.label == "SIGF":
            col.label = "F_SIGF_SIGF"
    src = tmp_path / "complete.mtz"
    mtz.write_to_file(str(src))
    root = tmp_path / "staging"
    manifest = stage_datasets([DatasetSpec(label="i2", xyzin=spec.xyzin, hklin=src)], root)
    staged = gemmi.read_mtz_file(str(root / DATASETS_DIR / xtal_name(0) / REFLECTIONS_NAME))
    labels = [c.label for c in staged.columns]
    assert "FreeR_flag" in labels and "FREERFLAG_FREER" not in labels
    assert [c.label for c in staged.columns if c.type == "I"] == ["FreeR_flag"]
    assert manifest["datasets"][0]["files"][REFLECTIONS_NAME]["how"] == "rewritten"


def test_dictionary_is_respelled_where_the_source_needs_it(tmp_path, mini_specs):
    spec = mini_specs[0]
    pdbx = tmp_path / "ligand_pdbx.cif"
    pdbx.write_text(pdbx_spelling(open(spec.dictionary).read()))
    specs = [DatasetSpec(label="needs-type", xyzin=spec.xyzin, hklin=spec.hklin, dictionary=pdbx)]
    root = tmp_path / "staging"
    manifest = stage_datasets(specs, root)
    staged = root / DATASETS_DIR / xtal_name(0) / DICT_NAME
    block = next(b for b in gemmi.cif.read(str(staged)) if b.find_values("_chem_comp_bond.atom_id_1"))
    assert block.find_values("_chem_comp_bond.type")
    assert manifest["datasets"][0]["files"][DICT_NAME]["how"] == "rewritten"
    assert not list((root / DATASETS_DIR / xtal_name(0)).glob("*_pandda.cif"))


def test_refuses_duplicate_labels_missing_files_and_existing_tree(tmp_path, mini_specs):
    with pytest.raises(ValueError, match="duplicate"):
        stage_datasets([mini_specs[0], mini_specs[0]], tmp_path / "a")
    with pytest.raises(ValueError):
        stage_datasets([], tmp_path / "b")
    bad = DatasetSpec(label="ghost", xyzin=tmp_path / "nope.pdb", hklin=mini_specs[0].hklin)
    with pytest.raises(FileNotFoundError, match="ghost"):
        stage_datasets([bad], tmp_path / "c")
    assert not (tmp_path / "c").exists(), "checked before anything is written"
    stage_datasets(mini_specs[:1], tmp_path / "d")
    with pytest.raises(FileExistsError):
        stage_datasets(mini_specs[:1], tmp_path / "d")


def test_link_or_copy_never_symlinks_and_resolves_a_symlinked_source(tmp_path):
    real = tmp_path / "real.bin"
    real.write_bytes(b"x" * 1000)
    link = tmp_path / "link.bin"
    link.symlink_to(real)
    dst = tmp_path / "out" / "staged.bin"
    how = link_or_copy(link, dst)
    assert how in ("link", "copy")
    assert not dst.is_symlink()
    assert dst.read_bytes() == real.read_bytes()
    if how == "link":
        assert os.stat(dst).st_ino == os.stat(real).st_ino, "linked to the file, not the link"
    with pytest.raises(FileExistsError):
        link_or_copy(real, dst)


@pytest.mark.skipif(not REAL_TREE.is_dir(), reason="BAZ2B staged tree not mounted")
def test_reproduces_the_real_staged_tree(tmp_path, mini_specs):
    """The first three datasets of the real 201-dataset tree, reproduced from
    their sources: same names, same Projects.csv rows, same bytes."""
    root = tmp_path / "staging"
    stage_datasets(mini_specs, root)
    real_rows = read_projects_csv(REAL_TREE / PROJECTS_CSV)[:3]
    assert read_projects_csv(root / PROJECTS_CSV) == real_rows
    for xtal, _label in real_rows:
        for name in (MODEL_NAME, REFLECTIONS_NAME, DICT_NAME):
            assert filecmp.cmp(REAL_TREE / DATASETS_DIR / xtal / name,
                               root / DATASETS_DIR / xtal / name, shallow=False), (xtal, name)
