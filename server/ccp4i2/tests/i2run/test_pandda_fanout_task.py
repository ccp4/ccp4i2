"""The fan-out task (design note section 8, decision 1 as amended): a job
that takes the manifest as a typed input and records, per dataset, what it
created. Receipts are created pending here (the test database is in memory,
so a subprocess could not see them); running them is covered elsewhere."""
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from django.conf import settings

from ccp4i2.db import models
from ccp4i2.tests.unit.pandda.synthetic_tree import event_record, make_tree
from .utils import i2run

pytest.importorskip("yaml", reason="needs PyYAML")


def _member(name):
    directory = Path(settings.CCP4I2_PROJECTS_DIR) / name
    (directory / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    return models.Project.objects.create(name=name, directory=str(directory))


def _fixture(tmp_path, members):
    tree = make_tree(tmp_path / "pandda2_out", {"xtal-0000": [event_record(1)], "xtal-0001": []})
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "manifest_version": 1, "created": "now", "datasets_dir": "datasets", "projects_csv": "Projects.csv",
        "provenance": {},
        "datasets": [
            {"xtal": "xtal-0000", "label": members[0].name, "project_uuid": str(members[0].uuid), "files": {}},
            {"xtal": "xtal-0001", "label": members[1].name, "project_uuid": str(members[1].uuid), "files": {}},
            {"xtal": "xtal-0002", "label": "absent", "project_uuid": str(members[1].uuid), "files": {}},
        ],
    }))
    return tree, manifest


def test_fanout_task_creates_receipts_in_the_members_projects(tmp_path):
    members = [_member("fan_m1"), _member("fan_m2")]
    tree, manifest = _fixture(tmp_path, members)
    args = ["pandda_fanout", "--MANIFEST", str(manifest), "--PANDDA_OUT_DIR", str(tree),
            "--RUN_RECEIPTS", "False"]
    with i2run(args, allow_errors=True) as job:
        job_id = ET.parse(job / "params.xml").find(".//jobId").text
        db_job = models.Job.objects.get(uuid=job_id)
        assert db_job.status == models.Job.Status.FINISHED
        receipts = models.Job.objects.filter(task_name="pandda_events").select_related("project")
        assert sorted(r.project.name for r in receipts) == ["fan_m1", "fan_m2"]
        assert all(r.parent_id is None for r in receipts), "receipts are top-level jobs of their own projects"
        xml = ET.parse(job / "program.xml")
        outcomes = {d.get("xtal"): d.get("action") for d in xml.findall("datasets/dataset")}
        assert outcomes == {"xtal-0000": "created", "xtal-0001": "created", "xtal-0002": "absent"}
        kpis = {v.key.name: v.value for v in models.JobFloatValue.objects.select_related("key").filter(job=db_job)}
        assert (kpis["nCreated"], kpis["nAbsent"]) == (2, 1)
        codes = {r.findtext("code") for r in ET.parse(job / "diagnostic.xml").findall(".//errorReport")}
        assert "204" in codes


def test_fanout_task_preview_creates_nothing(tmp_path):
    members = [_member("fan_p1"), _member("fan_p2")]
    tree, manifest = _fixture(tmp_path, members)
    args = ["pandda_fanout", "--MANIFEST", str(manifest), "--PANDDA_OUT_DIR", str(tree), "--DRY_RUN", "True"]
    with i2run(args, allow_errors=True) as job:
        assert models.Job.objects.filter(task_name="pandda_events").count() == 0
        xml = ET.parse(job / "program.xml")
        assert xml.findtext("dry_run") == "True" and xml.findtext("n_created") == "2"


def test_fanout_task_refuses_a_manifest_with_no_tree(tmp_path):
    members = [_member("fan_r1"), _member("fan_r2")]
    _tree, manifest = _fixture(tmp_path, members)
    with i2run(["pandda_fanout", "--MANIFEST", str(manifest)], allow_errors=True) as job:
        job_id = ET.parse(job / "params.xml").find(".//jobId").text
        assert models.Job.objects.get(uuid=job_id).status == models.Job.Status.FAILED
