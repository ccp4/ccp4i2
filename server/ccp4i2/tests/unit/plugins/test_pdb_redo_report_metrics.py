"""The PDB-REDO report shows what the run changed, before and after.

The metric tables were a commented-out block that never ran, so the report
showed links and logs and no numbers. The values below are from a real run
(2d2k, an RNA, 2026-10-07): the protein-only measures come back as None
there, and their rows must not appear.
"""

import xml.etree.ElementTree as ET

import pytest

django = pytest.importorskip("django")


@pytest.fixture(scope="module", autouse=True)
def _django_setup():
    import os
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.test_settings")
    django.setup()


def rendered_text(report):
    return " ".join(t.text or "" for t in report.as_data_etree().iter() if t.text)


RNA_RUN = {
    "PDB_REDO_JOB_ID": "2", "RCAL": "0.2472", "RFIN": "0.2362",
    "RFCAL": "0.2495", "RFFIN": "0.2477", "OBRMSZ": "0.835", "FBRMSZ": "0.463",
    "OARMSZ": "1.496", "FARMSZ": "0.814", "OCLASH": "32.76", "FCLASH": "11.87",
    "TOZRAMA": "None", "TFZRAMA": "None", "TOCONFAL": "68", "TFCONFAL": "72",
    "NBBFLIP": "0", "NWATDEL": "None", "RSCCB": "13", "RSCCW": "0",
}


def report_for(values):
    from ccp4i2.wrappers.pdb_redo_api.script.pdb_redo_api_report import pdb_redo_api_report

    root = ET.Element("pdb_redo_api")
    for tag, text in values.items():
        ET.SubElement(root, tag).text = text
    return pdb_redo_api_report(xmlnode=root, jobInfo={}, jobStatus=None)


def test_before_and_after_are_shown():
    text = rendered_text(report_for(RNA_RUN))
    assert "Refinement" in text and "Model changes" in text
    for value in ("0.2495", "0.2477", "32.76", "11.87", "68", "72", "13"):
        assert value in text, value
    assert "Dinucleotide conformation (CONFAL)" in text


def test_measures_the_run_did_not_report_are_left_out():
    text = rendered_text(report_for(RNA_RUN))
    assert "Ramachandran plot appearance" not in text
    assert "Waters deleted" not in text
    assert "None" not in text


def test_no_measures_no_tables():
    text = rendered_text(report_for({"PDB_REDO_JOB_ID": "2"}))
    assert "Refinement" not in text and "Model changes" not in text
