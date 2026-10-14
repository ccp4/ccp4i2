"""Looking at a run's events before the receipts exist: the list the report
offers, and the scene one of them opens. CCP4-free (gemmi for the maps)."""
import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

pytest.importorskip("gemmi", reason="needs gemmi")

import ccp4i2
from ccp4i2.wrappers.pandda_campaign.script.pandda_campaign_report import pandda_campaign_report
from ccp4i2.wrappers.pandda_campaign.script.pandda_live_scene import live_event_scene, peek_events

from .synthetic_tree import event_record, make_tree

CLIENT_SCHEMA = os.path.join(os.path.dirname(ccp4i2.__file__), "..", "..", "client", "renderer",
                             "lib", "scene", "moorhen-scene.ccp4i2.v1.json")


@pytest.fixture
def run(tmp_path):
    """A project directory holding a run: two datasets with events (one with
    two), one with none, one PanDDA has not finished."""
    tree = make_tree(tmp_path / "CCP4_JOBS" / "job_3" / "pandda2_out", {
        "xtal-0001": [event_record(1, score=0.2), event_record(2, score=0.9, build=False)],
        "xtal-0002": [event_record(1, score=0.5)],
        "xtal-0003": [],
    }, events_table=False)
    unfinished = tree / "processed_datasets" / "xtal-0004"
    unfinished.mkdir()
    (unfinished / "events.yaml").write_text("1: {Score: 0.99}\n")
    return tmp_path, tree


def test_peek_lists_finished_datasets_events_best_first(run):
    _, tree = run
    events = peek_events(tree)
    assert [(e["dtag"], e["event_idx"]) for e in events] == [
        ("xtal-0001", 2), ("xtal-0002", 1), ("xtal-0001", 1)]


def test_peek_reads_a_finished_dataset_once(run):
    _, tree = run
    cache = {}
    peek_events(tree, cache)
    assert set(cache) == {"xtal-0001", "xtal-0002", "xtal-0003"}
    (tree / "processed_datasets" / "xtal-0002" / "events.yaml").write_text("{}\n")
    assert len(peek_events(tree, cache)) == 3      # from the cache, not re-read


def test_event_scene_names_the_runs_files_by_project_path(run):
    project, tree = run
    scene = live_event_scene(tree, None, project, "uuid-1", "xtal-0001", 1)
    paths = {f["name"]: f.get("projectPath") for f in scene["files"]}
    assert paths["event1_map"] == ("CCP4_JOBS/job_3/pandda2_out/processed_datasets/xtal-0001/"
                                   "xtal-0001-event_1_1-BDC_0.2_map.native.ccp4")
    assert paths["zmap"].endswith("xtal-0001/xtal-0001-z_map.native.ccp4")
    assert all(f["projectId"] == "uuid-1" for f in scene["files"])
    maps = {m["file"]: m for m in scene["maps"]}
    assert maps["zmap"]["visible"] is False
    assert maps["event1_map"]["visible"] is True and maps["event1_map"]["contourLevel"] > 0
    assert scene["view"]["origin"] == [-15.0, -40.0, -30.0]


def test_files_outside_the_project_are_left_out(run):
    """The synthetic apo link points beside the tree (staged_inputs/), so
    with the tree as the project root it lies outside: the scene must not
    name it (the server would refuse to serve it anyway)."""
    _, tree = run
    scene = live_event_scene(tree, None, tree, "uuid-1", "xtal-0002", 1)
    assert "apo" not in {f["name"] for f in scene["files"]}


def test_scenes_validate_against_the_contract(run):
    jsonschema = pytest.importorskip("jsonschema", reason="needs jsonschema")
    if not os.path.isfile(CLIENT_SCHEMA):
        pytest.skip("client scene schema not beside this checkout")
    schema = json.loads(Path(CLIENT_SCHEMA).read_text())
    project, tree = run
    for dtag, idx in (("xtal-0001", 1), ("xtal-0001", 2), ("xtal-0003", None)):
        jsonschema.validate(live_event_scene(tree, None, project, "u", dtag, idx), schema)


@pytest.mark.parametrize("dtag,idx", [("../job_2", None), ("xtal-0001", 7), ("xtal-0009", 1)])
def test_bad_requests_say_why(run, dtag, idx):
    project, tree = run
    with pytest.raises(ValueError):
        live_event_scene(tree, None, project, "u", dtag, idx)


def test_the_report_offers_a_link_per_event():
    xml = ("<pandda_campaign><state>running</state>"
           "<datasets><dataset xtal='xtal-0001' label='BAZ2B-x1' dict='yes'/></datasets>"
           "<peek total='2' datasets='5'><event dtag='xtal-0001' idx='2' score='0.900'/>"
           "<event dtag='xtal-0002' idx='1'/></peek></pandda_campaign>")
    report = pandda_campaign_report(xmlnode=ET.fromstring(xml), jobInfo={"project_pk": 4},
                                    jobStatus="Running")
    links = [e for e in report.as_data_etree().iter() if e.tag == "CCP4i2ReportMoorhenLink"]
    assert [l.get("label") for l in links] == ["BAZ2B-x1 event 2 (score 0.900)", "xtal-0002 event 1"]
    assert json.loads(links[0].get("kwargs")) == {"dtag": "xtal-0001", "event_idx": 2}
    assert links[0].get("method") == "eventScene" and links[0].get("projectId") == "4"
