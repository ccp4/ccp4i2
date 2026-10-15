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
from ccp4i2.wrappers.pandda_campaign.script.pandda_dict import add_lig_alias
from ccp4i2.wrappers.pandda_campaign.script.pandda_live_scene import (
    dictionary_for_coot, live_event_scene, peek_events)

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
    assert all(f["projectId"] == "uuid-1" for f in scene["files"] if f["name"] != "dict")
    # The dictionary travels inline, rewritten for Coot (see below)
    dictionary = next(f for f in scene["files"] if f["name"] == "dict")
    assert "projectPath" not in dictionary and "data_comp_LIG" in dictionary["cifText"]
    maps = {m["file"]: m for m in scene["maps"]}
    assert maps["zmap"]["visible"] is False
    assert maps["event1_map"]["visible"] is True and maps["event1_map"]["contourLevel"] > 0
    assert scene["view"]["origin"] == [-15.0, -40.0, -30.0]
    dicts = {e["file"]: e.get("dictionaries") for e in scene["elements"]}
    assert dicts == {"apo": ["dict"], "event1_pose": ["dict"]}


def test_the_dictionary_rides_on_the_protein_without_a_pose(run):
    """Event 2 has no build: no pose, but the ligand's dictionary still goes
    on the apo model, for a ligand built there by hand."""
    project, tree = run
    scene = live_event_scene(tree, None, project, "u", "xtal-0001", 2)
    assert {f["name"] for f in scene["files"]} >= {"apo", "dict"}
    assert next(e for e in scene["elements"] if e["file"] == "apo")["dictionaries"] == ["dict"]


def test_the_apo_model_can_be_a_database_file(run):
    """On the campaign page the member's current model is already loaded:
    the scene names it by file id, so the viewer reuses it (and attaches
    the dictionary to it) instead of loading the staged apo model again."""
    project, tree = run
    scene = live_event_scene(tree, None, project, "uuid-1", "xtal-0001", 1, apo_file_id=321)
    apo = next(f for f in scene["files"] if f["name"] == "apo")
    assert apo == {"name": "apo", "kind": "coordinates", "fileId": 321, "projectId": "uuid-1"}
    assert next(e for e in scene["elements"] if e["file"] == "apo")["dictionaries"] == ["dict"]
    # the dataset-only view too
    scene = live_event_scene(tree, None, project, "uuid-1", "xtal-0003", None, apo_file_id=321)
    assert next(f for f in scene["files"] if f["name"] == "apo")["fileId"] == 321


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


def test_the_report_opens_a_campaign_members_event_on_the_campaign_page():
    """With the campaign and the member's current-model job recorded, the
    link opens the campaign viewer on that job; a member with no current
    model, or a run outside a campaign, opens on the project page."""
    xml = ("<pandda_campaign><state>finished</state>"
           "<peek total='2' datasets='5' campaign='9'>"
           "<event dtag='xtal-0001' idx='1' member_project='12' member_job='345'/>"
           "<event dtag='xtal-0002' idx='1' member_project='13'/></peek></pandda_campaign>")
    report = pandda_campaign_report(xmlnode=ET.fromstring(xml), jobInfo={"project_pk": 4},
                                    jobStatus="Finished")
    links = [e for e in report.as_data_etree().iter() if e.tag == "CCP4i2ReportMoorhenLink"]
    assert (links[0].get("campaignId"), links[0].get("jobId")) == ("9", "345")
    assert links[1].get("campaignId") is None and links[1].get("jobId") is None
    text = " ".join(e.text or "" for e in report.as_data_etree().iter() if (e.text or "").strip())
    assert "campaign viewer" in text


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


def test_the_dictionary_is_rewritten_so_coot_knows_lig(tmp_path):
    """Staging's comp_LIG alias is a copy whose rows still name the true
    code: PanDDA opens it by block name, Coot registers by comp_id and
    comp_list, so to Coot the pose's LIG had no dictionary and Moorhen
    fetched the library's unrelated LIG."""
    gemmi = pytest.importorskip("gemmi")
    from .synthetic_tree import DICT_CIF
    doc = gemmi.cif.read_string(DICT_CIF.format(code="MZ0"))
    add_lig_alias(doc)
    path = tmp_path / "dict.cif"
    path.write_text(doc.as_string())
    assert set(gemmi.cif.read(str(path))["comp_LIG"].find_values("_chem_comp_atom.comp_id")) == {"MZ0"}

    rewritten = gemmi.cif.read_string(dictionary_for_coot(path))
    assert list(rewritten["comp_list"].find_values("_chem_comp.id")) == ["MZ0", "LIG"]
    assert set(rewritten["comp_LIG"].find_values("_chem_comp_atom.comp_id")) == {"LIG"}
    assert set(rewritten["comp_MZ0"].find_values("_chem_comp_atom.comp_id")) == {"MZ0"}
    # the file itself is untouched
    assert set(gemmi.cif.read(str(path))["comp_LIG"].find_values("_chem_comp_atom.comp_id")) == {"MZ0"}


def test_a_dictionary_without_an_alias_gets_one(tmp_path):
    gemmi = pytest.importorskip("gemmi")
    from .synthetic_tree import DICT_CIF
    path = tmp_path / "dict.cif"
    path.write_text(DICT_CIF.format(code="MZ0"))
    rewritten = gemmi.cif.read_string(dictionary_for_coot(path))
    assert set(rewritten["comp_LIG"].find_values("_chem_comp_atom.comp_id")) == {"LIG"}
    assert "LIG" in list(rewritten["comp_list"].find_values("_chem_comp.id"))
