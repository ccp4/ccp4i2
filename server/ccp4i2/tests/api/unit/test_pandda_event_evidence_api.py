"""The evidence for one PanDDA event, as the campaign viewer asks for it.

    POST /api/ccp4i2/jobs/{receipt}/object_method/
         {"object_path": "pandda_events", "method_name": "eventEvidence",
          "args": [event_idx]}

A filled box in the campaign's site matrix is an event; clicking it opens the
dataset's model at the site and then asks the event's receipt for its event
map, autobuilt pose and dictionary through this plugin method. The tests hold
the two things that are easy to get wrong: an event is found by its PanDDA
number, not by its position in the receipt's list (``EVENTS[k]``), and the
contour is the display contour, falling back to the optimal one, in absolute
map units.

Uses the api/ conftest, which auto-applies django_db(transaction=True) and
sets AllowAny on the viewsets -- do NOT add @pytest.mark.django_db here.
"""
import uuid

from rest_framework.test import APIClient

from ccp4i2.db import models

from .test_pandda_site_index_api import _receipt


def _project(root, name="frag_ev"):
    root.mkdir(parents=True, exist_ok=True)
    return models.Project.objects.create(name=name, directory=str(root / name))


def _file(job, name, type_name, param_name):
    file_type, _ = models.FileType.objects.get_or_create(name=type_name)
    (job.directory / name).write_text("x")
    return models.File.objects.create(
        uuid=uuid.uuid4(), name=name, directory=models.File.Directory.JOB_DIR,
        type=file_type, job=job, job_param_name=param_name)


def _annotate(job, event_idx, extra):
    """Add fields the shared fixture does not write to one event's record."""
    params = job.directory / "params.xml"
    marker = f"<EVENT_IDX>{event_idx}</EVENT_IDX>"
    text = params.read_text()
    assert marker in text
    params.write_text(text.replace(marker, marker + extra))


def _receipt_with_files(root):
    """A receipt whose events are numbered 3 and 1, in that order: event 1
    sits at position 1 and event 3 at position 0, so a lookup that confused
    the two would fetch the other event's files."""
    project = _project(root)
    job = _receipt(project, "x0101", [
        {"idx": 3, "site": 2, "score": 0.2, "centroid": (1.0, 2.0, 3.0)},
        {"idx": 1, "site": 1, "score": 0.9, "centroid": (5.0, 5.0, 5.0)},
    ])
    _annotate(job, 1, "<OPTIMAL_CONTOUR>0.61</OPTIMAL_CONTOUR>"
                      "<DISPLAY_CONTOUR>0.42</DISPLAY_CONTOUR>"
                      "<LIGAND_ID>ZZY</LIGAND_ID>")
    _annotate(job, 3, "<OPTIMAL_CONTOUR>0.75</OPTIMAL_CONTOUR>")
    files = {
        "map0": _file(job, "event_3_map.map", "application/CCP4-map", "EVENTS[0].EVENT_MAP"),
        "pose0": _file(job, "event_3_pose.pdb", "chemical/x-pdb", "EVENTS[0].POSE"),
        "map1": _file(job, "event_1_map.map", "application/CCP4-map", "EVENTS[1].EVENT_MAP"),
        "pose1": _file(job, "event_1_pose.pdb", "chemical/x-pdb", "EVENTS[1].POSE"),
        "dict": _file(job, "ZZY.cif", "application/refmac-dictionary", "DICT"),
    }
    return job, files


def _evidence(job, event_idx):
    response = APIClient().post(
        f"/api/ccp4i2/jobs/{job.id}/object_method/",
        {"object_path": "pandda_events", "method_name": "eventEvidence",
         "args": [event_idx], "kwargs": {}},
        format="json")
    assert response.status_code == 200, response.content
    return response.json()["data"]["result"]


def test_an_event_is_found_by_its_number_not_its_position(
        bypass_api_permissions, test_project_path):
    job, files = _receipt_with_files(test_project_path)

    result = _evidence(job, 1)
    assert result["success"] is True, result
    data = result["data"]
    assert data["event_idx"] == 1 and data["position"] == 1
    assert data["receipt_job_id"] == job.id
    assert data["event_map"]["id"] == files["map1"].id
    assert data["event_map"]["uuid"] == str(files["map1"].uuid)
    assert data["event_map"]["type"] == "application/CCP4-map"
    assert data["pose"]["id"] == files["pose1"].id
    assert data["dictionary"]["id"] == files["dict"].id
    assert data["has_map"] is True and data["has_pose"] is True
    # The display contour wins over the optimal one.
    assert data["contour"] == 0.42
    assert data["display_contour"] == 0.42 and data["optimal_contour"] == 0.61
    assert data["centroid"] == [5.0, 5.0, 5.0]
    assert data["ligand_id"] == "ZZY"
    assert data["colour"] == "#3f51b5"


def test_the_optimal_contour_stands_in_when_there_is_no_display_contour(
        bypass_api_permissions, test_project_path):
    job, files = _receipt_with_files(test_project_path)
    data = _evidence(job, 3)["data"]
    assert data["position"] == 0
    assert data["event_map"]["id"] == files["map0"].id
    assert data["contour"] == 0.75 and data["display_contour"] is None
    assert data["ligand_id"] is None


def test_a_missing_pose_is_absent_not_an_error(
        bypass_api_permissions, test_project_path):
    job, files = _receipt_with_files(test_project_path)
    files["pose1"].delete()
    data = _evidence(job, 1)["data"]
    assert data["pose"] is None and data["has_pose"] is False
    assert data["has_map"] is True


def test_an_unknown_event_says_so(bypass_api_permissions, test_project_path):
    job, _ = _receipt_with_files(test_project_path)
    result = _evidence(job, 9)
    assert result["success"] is False
    assert "no event 9" in result["error"]
