"""Unit tests for the bibliography util (lib/utils/jobs/bibliography.py).

Covers MedLine parsing, task->citation-key expansion, dedup, and that every
citation key referenced by the map / task identities resolves to a real file
(or is explicitly non-citable). No CCP4 binaries or DB needed.
"""

import pytest

from ccp4i2.lib.utils.jobs import bibliography as bib


def test_parse_medline_basic():
    text = (
        "\nPMID- 19461840\n"
        "TI  - Phaser crystallographic software.\n"
        "SO  - J Appl Crystallogr. 2007;40(Pt 4):658-674.\n"
        "AU  - McCoy AJ\n"
        "AU  - Read RJ\n"
        "URL - https://doi.org/10.1107/S0021889807021206\n"
    )
    refs = bib._parse_medline(text)
    assert len(refs) == 1
    r = refs[0]
    assert r["pmid"] == "19461840"
    assert r["title"] == "Phaser crystallographic software."
    assert r["authors"] == ["McCoy AJ", "Read RJ"]
    assert r["source"].startswith("J Appl Crystallogr")
    assert r["link"] == "https://doi.org/10.1107/S0021889807021206"


def test_parse_medline_doi_fallback_link():
    # No URL/UR line: link should be synthesised from the AID/LID doi.
    text = (
        "\nPMID- 21460454\n"
        "TI  - REFMAC5 for the refinement of macromolecular crystal structures.\n"
        "SO  - Acta Crystallogr D Biol Crystallogr. 2011;67(Pt 4):355-67.\n"
        "AU  - Murshudov GN\n"
        "AID - 10.1107/S0907444911001314 [doi]\n"
    )
    refs = bib._parse_medline(text)
    assert refs[0]["link"] == "https://doi.org/10.1107/S0907444911001314"


def test_parse_medline_multiline_title():
    text = (
        "\nPMID- 1\n"
        "TI  - A very long title that wraps across\n"
        "      a second indented line.\n"
        "SO  - Journal. 2020.\n"
    )
    refs = bib._parse_medline(text)
    assert refs[0]["title"] == "A very long title that wraps across a second indented line."


def test_citation_key_identity_and_alias():
    # Unmapped task cites itself.
    assert bib._citation_keys_for_task("refmac") == ["refmac"]
    # Variant aliases to canonical program.
    assert bib._citation_keys_for_task("phaser_MR_FRF") == ["phaser"]
    # Monolithic pipeline expands to its internal toolchain.
    assert set(bib._citation_keys_for_task("modelcraft")) >= {
        "modelcraft", "buccaneer", "nautilus", "sheetbend", "parrot", "refmac",
    }


def test_references_for_tasks_dedup():
    # crank2 pulls many keys; the CCP4-suite overview appears once even if
    # multiple keys cite it, and every ref has a dedup-distinct identity.
    refs = bib.references_for_tasks(["crank2"])
    assert len(refs) > 5
    keys = [bib._dedup_key(r) for r in refs]
    assert len(keys) == len(set(keys)), "references_for_tasks returned duplicates"


def test_refmac_has_expected_reference():
    refs = bib.references_for_tasks(["refmac"])
    titles = " ".join((r["title"] or "") for r in refs)
    assert "REFMAC5" in titles


def test_non_citable_returns_empty():
    assert bib.references_for_tasks(["ImportObs"]) == []
    assert bib.references_for_tasks(["splitMtz"]) == []


def test_every_cited_key_resolves():
    """Every citation key produced by the map (or a task identity) must resolve
    to a real MedLine file, unless the key is explicitly non-citable. Guards
    against a typo'd alias silently yielding no reference."""
    from ccp4i2.core.tasks import TASKS

    missing = []
    for task in TASKS:
        for key in bib._citation_keys_for_task(task):
            path = bib.REFERENCES_DIR / f"{key}.medline.txt"
            if not path.exists() and key not in bib._NON_CITABLE:
                missing.append((task, key))
    assert not missing, f"citation keys with no reference file: {missing}"


def test_all_reference_files_parse():
    """Every *.medline.txt in the references dir parses to >=1 reference (no
    silently-empty file from a bad edit)."""
    empty = []
    for path in bib.REFERENCES_DIR.glob("*.medline.txt"):
        text = path.read_text(encoding="utf-8", errors="replace")
        if not bib._parse_medline(text):
            empty.append(path.name)
    assert not empty, f"reference files that parse to zero references: {empty}"
