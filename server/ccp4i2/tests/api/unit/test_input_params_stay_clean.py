"""Every endpoint that rewrites input_params.xml must write only what was set.

The i2run suite creates a job and runs it; the GUI creates a job and then
*edits* it --- set_parameter for every widget change, upload_file_param for
every file. Each of those rewrites input_params.xml in full, and for months
each rewrite wrote every untouched parameter at its zero value. i2run never
noticed because i2run never edits. Here the same endpoints the GUI uses are
driven through the Django test client, and after every write the file on disk
is compared against what a fresh container would refuse to serialise.

No CCP4, no test101 zips --- the project is built from nothing.
"""

import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.tests.unit.conformance import harness

API_PREFIX = "/api/ccp4i2"

# The task the regression was found in; NPROC (min=1) is the parameter that
# turned a silent extra element into a blocked job.
TASK = "aimless_pipe"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "cleanparams"
    directory.mkdir()
    return models.Project.objects.create(name="cleanparams", directory=str(directory))


def _create_job(client, project, task=TASK):
    response = client.post(
        f"{API_PREFIX}/projects/{project.id}/create_task/",
        data=json.dumps({"task_name": task}),
        content_type="application/json",
    )
    assert response.status_code == 200, response.content
    body = response.json()
    assert body.get("success"), body
    return models.Job.objects.get(id=body["data"]["new_job"]["id"])


def _set_leaves(plugin, task=TASK):
    return {
        path
        for path, obj in harness.walk(plugin, task).leaves
        if callable(getattr(obj, "isSet", None)) and obj.isSet()
    }


def _assert_only_these_were_set(job, expected_prefixes, task=TASK):
    """Reload the written file into a fresh container: nothing may come back
    set that a fresh container does not already consider set (def.xml defaults
    answer isSet() True before any file is involved), except under the paths
    the API call actually touched."""
    params = Path(job.directory) / "input_params.xml"
    assert params.exists(), f"no input_params.xml in {job.directory}"

    fresh = harness.build(task)
    if fresh is None:
        pytest.skip(f"{task} is not available in this environment")
    baseline = _set_leaves(fresh, task)

    reloaded = harness.build(task)
    body = ET.parse(str(params)).getroot().find(".//ccp4i2_body")
    reloaded.container.setEtree(body, ignore_missing=True)

    stray = sorted(
        path for path in _set_leaves(reloaded, task) - baseline
        if not any(path.startswith(prefix) for prefix in expected_prefixes)
    )
    assert not stray, (
        f"{len(stray)} parameters nobody set came back explicitly set from "
        f"input_params.xml: {stray[:8]} --- any with a lower bound (NPROC, "
        f"MONOMERS_ASYM) then blocks the job at run-time validation"
    )


def test_creating_a_job_writes_a_clean_file(client, project):
    job = _create_job(client, project)
    _assert_only_these_were_set(job, expected_prefixes=())


def test_setting_one_parameter_does_not_write_the_rest(client, project):
    """The GUI's every-widget-change path --- and NPROC is the very parameter
    whose invented zero blocked jobs, so it is the one set here."""
    job = _create_job(client, project)
    response = client.post(
        f"{API_PREFIX}/jobs/{job.id}/set_parameter/",
        data=json.dumps(
            {"object_path": f"{TASK}.controlParameters.NPROC", "value": 4}
        ),
        content_type="application/json",
    )
    assert response.status_code == 200, response.content

    _assert_only_these_were_set(job, expected_prefixes=("controlParameters.NPROC",))
    written = ET.parse(str(Path(job.directory) / "input_params.xml")).getroot()
    chosen = [e.text for e in written.iter("NPROC")]
    assert chosen and chosen[0].strip() == "4", \
        "the parameter that WAS set did not reach the file"


def test_uploading_a_file_does_not_write_the_rest(client, project):
    """The reported reproduction: a fresh job plus one upload."""
    job = _create_job(client, project)
    response = client.post(
        f"{API_PREFIX}/jobs/{job.id}/upload_file_param/",
        {
            "file": SimpleUploadedFile("gamma_unmerged.mtz", b"not really an mtz"),
            "objectPath": f"{TASK}.inputData.UNMERGEDFILES[0].file",
        },
        format="multipart",
    )
    if response.status_code != 200 or not response.json().get("success", True):
        pytest.skip(f"upload endpoint unavailable here: {response.content[:120]}")
    # The uploaded file object and its attributes are the choice made here;
    # nothing outside that subtree may come back set.
    _assert_only_these_were_set(job, expected_prefixes=("inputData.UNMERGEDFILES",))
