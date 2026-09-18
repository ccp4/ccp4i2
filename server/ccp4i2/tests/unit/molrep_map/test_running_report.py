"""molrep_map reports on the first hand while the second runs: the wrapper
writes its XML after each hand with the rest marked pending, and the report
renders that honestly."""

import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pytest

django = pytest.importorskip("django")


@pytest.fixture(scope="module", autouse=True)
def _django_setup():
    import os
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.test_settings")
    django.setup()


def rendered_text(report):
    return " ".join(t.text or "" for t in report.as_data_etree().iter() if t.text)


def test_wrapper_writes_a_partial_xml_after_the_first_hand(tmp_path):
    from ccp4i2.wrappers.molrep_map.script.molrep_map import molrep_map

    plugin = molrep_map(parent=None, workDirectory=str(tmp_path))
    plugin._results = {"Original": SimpleNamespace(placed=True, timed_out=False, score=0.0412, doc_path=None)}
    plugin._cc = {"Original": 0.61}
    plugin._write_program_xml(pending=["Inverted"])

    root = ET.parse(str(tmp_path / "program.xml")).getroot()
    assert root.get("pending") == "Inverted"
    assert root.find("recommendation") is None, "no verdict until both hands are done"
    assert root.find("Original").get("placed") == "true"
    assert root.find("Original").get("map_model_cc") == "0.6100"
    assert root.find("Inverted").get("pending") == "true"
    assert root.find("Inverted").get("placed") is None


def test_wrapper_final_xml_has_no_pending_marker(tmp_path):
    from ccp4i2.wrappers.molrep_map.script.molrep_map import molrep_map

    plugin = molrep_map(parent=None, workDirectory=str(tmp_path))
    plugin._results = {h: SimpleNamespace(placed=True, timed_out=False, score=0.04, doc_path=None)
                       for h in ("Original", "Inverted")}
    plugin._cc = {"Original": 0.61, "Inverted": 0.32}
    plugin._recommended, plugin._confidence = "Original", "confident"
    plugin._write_program_xml()
    root = ET.parse(str(tmp_path / "program.xml")).getroot()
    assert root.get("pending") is None
    assert root.find("recommendation").get("hand") == "Original"
    assert root.find("Inverted").get("placed") == "true"


def test_report_renders_the_running_state():
    from ccp4i2.wrappers.molrep_map.script.molrep_map_report import molrep_map_report

    assert molrep_map_report.RUNNING is True
    xml = ET.fromstring(
        '<molrep_map pending="Inverted">'
        '<Original placed="true" timed_out="false" score="0.0412" map_model_cc="0.6100"/>'
        '<Inverted pending="true"/>'
        '</molrep_map>')
    text = rendered_text(molrep_map_report(xmlnode=xml, jobInfo={}, jobStatus="Running"))
    assert "original hand is placed" in text
    assert "inverted hand is still running" in text
    assert "recommended" not in text


def test_report_renders_the_finished_state_as_before():
    from ccp4i2.wrappers.molrep_map.script.molrep_map_report import molrep_map_report

    xml = ET.fromstring(
        '<molrep_map><recommendation hand="Original" confidence="confident" cc_original="0.6100" cc_inverted="0.3200"/>'
        '<Original placed="true" timed_out="false" score="0.0412"/>'
        '<Inverted placed="true" timed_out="false" score="0.0398"/>'
        '</molrep_map>')
    text = rendered_text(molrep_map_report(xmlnode=xml, jobInfo={}, jobStatus="Finished"))
    assert "still running" not in text
    assert "clear winner" in text
