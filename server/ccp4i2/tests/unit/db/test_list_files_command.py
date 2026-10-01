"""list_files --job takes a job number within --project, and lists what the
job made as well as what it used.

It refused a number ("Job must be specified by UUID") even with the project
given, which its own help invited, and listed only the files the job used.
The JSON says what each file is to the job: parameter and annotation.
"""
import json
import uuid
from io import StringIO

import pytest
from django.core.management import call_command

from ccp4i2.db import models

pytestmark = pytest.mark.django_db


def _project(tmp_path):
    project = models.Project.objects.create(name="lf", directory=str(tmp_path))
    kind, _ = models.FileType.objects.get_or_create(name="chemical/x-pdb")
    made_by = models.Job.objects.create(project=project, number="1", title="refine")
    user = models.Job.objects.create(project=project, number="2", title="deposit")
    model = models.File.objects.create(
        uuid=uuid.uuid4(), name="XYZOUT.pdb", type=kind, job=made_by,
        job_param_name="XYZOUT", annotation="Model from refinement",
        directory=models.File.Directory.JOB_DIR)
    models.File.objects.create(
        uuid=uuid.uuid4(), name="MMCIFOUT.cif", type=kind, job=user,
        job_param_name="MMCIFOUT", annotation="Coordinates for upload",
        directory=models.File.Directory.JOB_DIR)
    models.FileUse.objects.create(file=model, job=user, job_param_name="XYZIN",
                                  role=models.FileUse.Role.IN)
    return project


def _files(*args):
    out = StringIO()
    call_command("list_files", *args, "--format", "json", stdout=out)
    return json.loads(out.getvalue())


def test_a_job_by_number_lists_what_it_made_and_used(tmp_path):
    _project(tmp_path)
    files = _files("--project", "lf", "--job", "2")
    assert {(f["job"]["number"], f["param"], f["annotation"]) for f in files} == {
        ("1", "XYZOUT", "Model from refinement"),
        ("2", "MMCIFOUT", "Coordinates for upload"),
    }


def test_a_number_without_a_project_says_what_is_needed(tmp_path):
    _project(tmp_path)
    out = StringIO()
    call_command("list_files", "--job", "2", stdout=out)
    assert "with --project" in out.getvalue()
