"""The receipt's scenes validate against the scene contract the server
carries (a generated mirror of the client's Zod schema). CCP4-free."""
import json
import os

import pytest

import ccp4i2
from ccp4i2.wrappers.pandda_events.script import pandda_scene

jsonschema = pytest.importorskip("jsonschema", reason="needs jsonschema")

# The ccp4i2-dialect schema (job+param references) lives in the client tree;
# the server mirror carries only the structured (LLM-output) shape, in which
# every key is required. In a repo checkout the client tree is beside us.
CLIENT_SCHEMA = os.path.join(os.path.dirname(ccp4i2.__file__), "..", "..", "client", "renderer",
                             "lib", "scene", "moorhen-scene.ccp4i2.v1.json")


def _schema():
    if not os.path.isfile(CLIENT_SCHEMA):
        pytest.skip("client scene schema not beside this checkout")
    with open(CLIENT_SCHEMA) as handle:
        return json.load(handle)


EVENT = {"idx": 1, "position": 0, "centroid": (15.0, 40.0, 30.0), "display_contour": 0.42,
         "score": 0.3, "ligand_id": "MZ0", "has_map": True, "has_pose": True}
EVENT2 = {"idx": 2, "position": 1, "centroid": (1.0, 2.0, 3.0), "display_contour": None,
          "score": 0.8, "ligand_id": "MZ0", "has_map": True, "has_pose": False}


def test_event_and_overview_scenes_validate():
    schema = _schema()
    kw = dict(apo=True, zmap=True, dictionary=True)
    for scene in (pandda_scene.event_scene("xtal-0004", "2", "uuid", EVENT, **kw),
                  pandda_scene.overview_scene("xtal-0004", "2", "uuid", [EVENT, EVENT2], **kw)):
        jsonschema.validate(scene, schema)


def test_scene_says_what_the_receipt_knows():
    scene = pandda_scene.event_scene("xtal-0004", "2", "uuid", EVENT, apo=True, zmap=True, dictionary=True)
    assert scene["scene"] == "xtal-0004 event 1 (MZ0)"
    maps = {m["name"]: m for m in scene["maps"]}
    assert maps["Z-map"]["contourLevel"] == 3.0 and maps["Z-map"]["isDifference"] is True
    assert maps["Event 1 map"]["contourLevel"] == 0.42
    assert scene["activeMap"] == "Event 1 map"
    assert scene["view"] == {"origin": [15.0, 40.0, 30.0], "zoom": pandda_scene.EVENT_ZOOM}
    pose = next(e for e in scene["elements"] if e["file"] == "event1_pose")
    assert pose["dictionaries"] == ["dict"]
    assert all(f["job"] == 2 and f["projectId"] == "uuid" for f in scene["files"])


def test_nothing_is_referenced_that_the_receipt_did_not_get():
    scene = pandda_scene.overview_scene("x", "3", "u", [EVENT2], apo=False, zmap=False, dictionary=False)
    assert [f["name"] for f in scene["files"]] == ["event2_map"]
    assert scene["elements"] == []
    assert "contourLevel" not in scene["maps"][0]
    assert pandda_scene.dump(scene).startswith("scene:")


def test_the_overview_opens_on_the_best_scored_event():
    scene = pandda_scene.overview_scene("x", "3", "u", [EVENT, EVENT2], apo=True, zmap=True, dictionary=True)
    maps = {m["name"]: m for m in scene["maps"]}
    assert maps["Event 2 map"]["visible"] is True and maps["Event 1 map"]["visible"] is False
    assert scene["view"]["origin"] == [1.0, 2.0, 3.0] and scene["activeMap"] == "Event 2 map"
    assert pandda_scene.focus_event([]) is None
    unscored = [dict(EVENT, score=None), dict(EVENT2, score=None)]
    assert pandda_scene.focus_event(unscored)["idx"] == 1
