"""``GET /api/ccp4i2/jobs/<id>/directory/``: one job's own tree.

The Directory tab used to fetch the whole project and walk down to the job.
On a campaign parent whose PanDDA job holds a directory per dataset that was
thousands of entries per poll to render a few.
"""
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def jobs(test_project_path):
    test_project_path.mkdir(exist_ok=True)
    project = models.Project.objects.create(
        name="jobdir", directory=str(test_project_path / "jobdir"))
    one = models.Job.objects.create(project=project, number="1", task_name="pandda_campaign",
                                    title="p", status=models.Job.Status.RUNNING_REMOTELY)
    two = models.Job.objects.create(project=project, number="2", task_name="i2Dimple",
                                    title="d", status=models.Job.Status.FINISHED)
    one.directory.mkdir(parents=True)
    (one.directory / "program.xml").write_text("<x/>")
    (one.directory / "pandda2_out" / "processed_datasets" / "xtal-0000").mkdir(parents=True)
    (one.directory / "pandda2_out" / "processed_datasets" / "xtal-0000" / "z_map.ccp4").write_bytes(b"m")
    two.directory.mkdir(parents=True)
    (two.directory / "final.mtz").write_bytes(b"MTZ ")
    return {"project": project, "one": one, "two": two}


def _names(container):
    return {node["name"] for node in container}


def test_lists_only_that_job(client, jobs):
    resp = client.get(f"/api/ccp4i2/jobs/{jobs['one'].id}/directory/")
    assert resp.status_code == 200, resp.content
    body = resp.json()
    assert body["status"] == "Success"
    assert _names(body["container"]) == {"program.xml", "pandda2_out"}

    resp = client.get(f"/api/ccp4i2/jobs/{jobs['two'].id}/directory/")
    assert _names(resp.json()["container"]) == {"final.mtz"}


def test_descends_into_the_job_tree(client, jobs):
    body = client.get(f"/api/ccp4i2/jobs/{jobs['one'].id}/directory/").json()
    out = next(n for n in body["container"] if n["name"] == "pandda2_out")
    processed = next(n for n in out["contents"] if n["name"] == "processed_datasets")
    assert _names(processed["contents"]) == {"xtal-0000"}


def test_a_job_with_no_directory_yet_lists_empty(client, jobs):
    job = models.Job.objects.create(project=jobs["project"], number="3", task_name="i2Dimple",
                                    title="new", status=models.Job.Status.PENDING)
    body = client.get(f"/api/ccp4i2/jobs/{job.id}/directory/").json()
    assert body["container"] == []


def test_unknown_job_is_404(client, jobs):
    assert client.get("/api/ccp4i2/jobs/999999/directory/").status_code == 404


# ?path= lists one folder: the lazy browser opens folders one at a time
# instead of fetching (and stat'ing) the whole tree on every poll.

def test_path_lists_one_level(client, jobs):
    body = client.get(f"/api/ccp4i2/jobs/{jobs['one'].id}/directory/?path=").json()
    assert _names(body["container"]) == {"program.xml", "pandda2_out"}
    out = next(n for n in body["container"] if n["name"] == "pandda2_out")
    assert out["type"] == "directory"
    assert "contents" not in out
    xml = next(n for n in body["container"] if n["name"] == "program.xml")
    assert xml["size"] == 4 and "mtime" in xml


def test_path_lists_a_subfolder(client, jobs):
    url = f"/api/ccp4i2/jobs/{jobs['one'].id}/directory/?path=pandda2_out/processed_datasets/xtal-0000"
    body = client.get(url).json()
    assert [n["name"] for n in body["container"]] == ["z_map.ccp4"]
    assert body["container"][0]["path"].endswith("xtal-0000/z_map.ccp4")


def test_path_may_not_leave_the_job(client, jobs):
    resp = client.get(f"/api/ccp4i2/jobs/{jobs['one'].id}/directory/?path=../job_2")
    assert resp.status_code == 400


def test_path_through_a_symlink_out_of_the_tree_is_refused(client, jobs, tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    (jobs["one"].directory / "escape").symlink_to(outside)
    resp = client.get(f"/api/ccp4i2/jobs/{jobs['one'].id}/directory/?path=escape")
    assert resp.status_code == 400


def test_path_of_a_job_with_no_directory_yet_lists_empty(client, jobs):
    job = models.Job.objects.create(project=jobs["project"], number="3", task_name="i2Dimple",
                                    title="new", status=models.Job.Status.PENDING)
    body = client.get(f"/api/ccp4i2/jobs/{job.id}/directory/?path=").json()
    assert body["container"] == []


def test_search_finds_names_anywhere_in_the_job(client, jobs):
    body = client.get(f"/api/ccp4i2/jobs/{jobs['one'].id}/directory/?search=Z_MAP").json()
    found = body["container"]
    assert found["truncated"] is False
    assert [i["relative"] for i in found["items"]] == [
        "pandda2_out/processed_datasets/xtal-0000/z_map.ccp4"]


def test_project_directory_lists_one_level(client, jobs):
    project = jobs["project"]
    body = client.get(f"/api/ccp4i2/projects/{project.id}/directory/?path=CCP4_JOBS").json()
    assert _names(body["container"]) == {"job_1", "job_2"}
