"""Fetching maps from EMDB: entry parsing, archive URLs, annotations and the
streamed gunzip download. CCP4-free; the network is replaced by fixtures
captured from the EMDB API on 2026-09-17."""

import gzip
import io
import json
import pathlib
import urllib.request

import pytest

from ccp4i2.lib.utils.files import repository_fetch as repo

FIXTURES = pathlib.Path(__file__).parent / "fixtures" / "emdb"


def entry(name):
    return json.loads((FIXTURES / f"{name}.json").read_text())


@pytest.mark.parametrize("text", ["11638", "EMD-11638", "emd-11638", "emd_11638", " EMD11638 "])
def test_entry_ids_normalise(text):
    assert repo.normalise_emdb_entry(text) == "EMD-11638"


@pytest.mark.parametrize("text", ["", "1cbs", "EMD-", "emd_abc"])
def test_bad_entry_ids_are_refused(text):
    with pytest.raises(repo.RepositoryError) as err:
        repo.normalise_emdb_entry(text)
    assert err.value.status == 400


def test_entry_with_everything_lists_main_halves_and_mask():
    files = repo.emdb_files(entry("EMD-11638"))
    assert [(f["kind"], f["file"], f["sub_type"], f["directory"]) for f in files] == [
        ("map", "emd_11638.map.gz", 1, "map"),
        ("half_map", "emd_11638_half_map_1.map.gz", 5, "other"),
        ("half_map", "emd_11638_half_map_2.map.gz", 5, "other"),
        ("mask", "emd_11638_msk_1.map", 4, "masks"),
    ]
    assert [f["label"] for f in files] == ["Main map", "Half map 1", "Half map 2", "Mask 1"]
    assert files[1]["pixel_spacing"] == pytest.approx(0.5332)
    assert files[1]["dimensions"] == [256, 256, 256]
    assert files[0]["contour_level"] == pytest.approx(0.116)


def test_entry_without_mask_or_without_halves_lists_fewer():
    kinds = lambda name: [f["kind"] for f in repo.emdb_files(entry(name))]
    assert kinds("EMD-8117") == ["map", "half_map", "half_map"]
    assert kinds("EMD-30210") == ["map"]


def test_summary_carries_title_resolution_and_fitted_models():
    summary = repo.emdb_entry_summary(entry("EMD-11638"))
    assert summary["entry"] == "EMD-11638"
    assert summary["resolution"] == pytest.approx(1.22)
    assert "apoferritin" in summary["title"].lower()
    assert summary["pdb_ids"] == ["7a4m"]
    assert len(summary["files"]) == 4


def test_archive_url_comes_from_the_kind_not_the_client():
    files = {f["file"]: f for f in repo.emdb_files(entry("EMD-11638"))}
    assert repo.emdb_archive_url("EMD-11638", files["emd_11638_half_map_2.map.gz"]) == \
        "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-11638/other/emd_11638_half_map_2.map.gz"
    assert repo.emdb_archive_url("EMD-11638", files["emd_11638_msk_1.map"]).endswith("/masks/emd_11638_msk_1.map")
    with pytest.raises(repo.RepositoryError):
        repo.emdb_archive_url("EMD-11638", {"file": "../../etc/passwd", "directory": "map"})


def test_annotation_says_what_and_how_fine():
    files = repo.emdb_files(entry("EMD-11638"))
    assert repo.emdb_annotation("EMD-11638", files[1], entry("EMD-11638")) == \
        "EMD-11638 half map 1, 0.53 A/px, 256^3, 1.22 A"
    assert repo.emdb_annotation("EMD-30210", repo.emdb_files(entry("EMD-30210"))[0], entry("EMD-30210")).startswith("EMD-30210 main map, 1.01 A/px, 192^3")


def test_download_streams_and_gunzips(tmp_path):
    payload = b"MAP " + bytes(range(256)) * 64
    source = tmp_path / "emd_1.map.gz"
    with gzip.open(source, "wb") as handle:
        handle.write(payload)
    target = tmp_path / "emd_1.map"
    repo.download_to(source.as_uri(), target, gunzip=True)
    assert target.read_bytes() == payload
    plain = tmp_path / "emd_2.map"
    plain.write_bytes(payload)
    assert repo.download_to(plain.as_uri(), tmp_path / "copy.map", gunzip=False).read_bytes() == payload


def test_download_failure_is_a_repository_error(tmp_path):
    with pytest.raises(repo.RepositoryError) as err:
        repo.download_to((tmp_path / "missing.map").as_uri(), tmp_path / "out.map", gunzip=False)
    assert err.value.status == 502
