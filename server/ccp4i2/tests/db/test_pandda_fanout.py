"""Fan-out (design note section 8, assertions in 15.2): twice over the same
manifest creates nothing the second time; over a partial tree it creates
receipts for what is there and reports the rest absent; against a tree with
no run job of its own it works, which is 5.4's second requirement as an
executable assertion. Receipts are created and parameterised through the
same helpers the interface uses; running them is the i2run tests' job."""

import json
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path
from shutil import rmtree

import pytest
from django.core.management import call_command
from django.test import TestCase, override_settings

from ...db import models
from ...lib.utils.jobs.pandda_fanout import plan_fanout
from ...tests.unit.pandda.synthetic_tree import event_record, make_tree

pytest.importorskip("yaml", reason="needs PyYAML")

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_FANOUT_TEST_DIR"


def make_project(name):
    directory = PROJECTS_DIR / name
    (directory / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    return models.Project.objects.create(name=name, directory=str(directory))


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class FanoutTest(TestCase):
    def setUp(self):
        PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)
        self.p0 = make_project("BAZ2BA-x425")
        self.p1 = make_project("BAZ2BA-x427")
        self.p2 = make_project("BAZ2BA-x428")
        self.tree = make_tree(PROJECTS_DIR / "run" / "pandda2_out", {
            "xtal-0000": [event_record(1)],
            "xtal-0001": [],
        })
        self.manifest = PROJECTS_DIR / "run" / "manifest.json"
        self.manifest.write_text(json.dumps({
            "manifest_version": 1, "created": "now", "datasets_dir": "datasets",
            "projects_csv": "Projects.csv",
            "provenance": {"run_job_uuid": "11111111-2222-3333-4444-555555555555"},
            "datasets": [
                {"xtal": "xtal-0000", "label": "BAZ2BA-x425", "project_uuid": str(self.p0.uuid), "files": {}},
                {"xtal": "xtal-0001", "label": "BAZ2BA-x427", "project_uuid": str(self.p1.uuid), "files": {}},
                {"xtal": "xtal-0002", "label": "BAZ2BA-x428", "project_uuid": str(self.p2.uuid), "files": {}},
                {"xtal": "xtal-0003", "label": "gone", "project_uuid": "00000000-0000-0000-0000-000000000000", "files": {}},
            ],
        }))

    def fanout(self, *extra):
        out = StringIO()
        call_command("pandda_fanout", "--tree", str(self.tree), "--manifest", str(self.manifest),
                     "--no-run", "--json", *extra, stdout=out)
        return json.loads(out.getvalue())

    @staticmethod
    def receipt_inputs(job):
        root = ET.parse(Path(job.directory) / "input_params.xml").getroot()
        inputs = root.find("ccp4i2_body/inputData")
        return {tag: (inputs.findtext(tag) or "") for tag in
                ("PANDDA_OUT_DIR", "DTAG", "RUN_JOB_UUID", "RUN_INCOMPLETE")}

    def test_a_dry_run_reports_the_plan_and_creates_nothing(self):
        report = self.fanout("--dry-run")
        self.assertTrue(report["dry_run"])
        self.assertEqual({o["xtal"]: o["action"] for o in report["outcomes"]}, {
            "xtal-0000": "created", "xtal-0001": "created",
            "xtal-0002": "absent", "xtal-0003": "no_project"})
        self.assertEqual(models.Job.objects.count(), 0)

    def test_receipts_land_in_their_own_projects_parameterised(self):
        report = self.fanout()
        actions = {o["xtal"]: o["action"] for o in report["outcomes"]}
        self.assertEqual(actions["xtal-0000"], "created")
        self.assertEqual(actions["xtal-0001"], "created")
        self.assertEqual(actions["xtal-0002"], "absent")
        self.assertEqual(models.Job.objects.filter(task_name="pandda_events").count(), 2)
        job = models.Job.objects.get(project=self.p0)
        self.assertEqual(job.status, models.Job.Status.PENDING)
        self.assertIn("BAZ2BA-x425", job.title)
        inputs = self.receipt_inputs(job)
        self.assertEqual(inputs["DTAG"], "xtal-0000")
        self.assertEqual(Path(inputs["PANDDA_OUT_DIR"]), self.tree)
        self.assertEqual(inputs["RUN_JOB_UUID"].replace("-", ""), "11111111222233334444555555555555")
        self.assertEqual(inputs["RUN_INCOMPLETE"], "False")
        self.assertEqual(models.Job.objects.filter(project=self.p2).count(), 0)

    def test_running_twice_creates_nothing_the_second_time(self):
        self.fanout()
        report = self.fanout()
        self.assertEqual({o["xtal"]: o["action"] for o in report["outcomes"]}, {
            "xtal-0000": "skipped", "xtal-0001": "skipped",
            "xtal-0002": "absent", "xtal-0003": "no_project"})
        self.assertEqual(models.Job.objects.filter(task_name="pandda_events").count(), 2)

    def test_a_partial_tree_gives_receipts_for_what_is_there_marked_incomplete(self):
        rmtree(self.tree / "analyses")
        report = self.fanout()
        self.assertTrue(report["incomplete"])
        self.assertEqual(report["outcomes"][0]["action"], "created")
        job = models.Job.objects.get(project=self.p0)
        self.assertEqual(self.receipt_inputs(job)["RUN_INCOMPLETE"], "True")

    def test_a_tree_produced_elsewhere_needs_no_run_job(self):
        data = json.loads(self.manifest.read_text())
        data["provenance"] = {}
        self.manifest.write_text(json.dumps(data))
        report = self.fanout()
        self.assertIsNone(report["run_job_uuid"])
        job = models.Job.objects.get(project=self.p0)
        self.assertEqual(self.receipt_inputs(job)["RUN_JOB_UUID"], "")
        self.assertIn("external", job.title)
        # still idempotent, keyed on the tree instead
        again = self.fanout()
        self.assertEqual(again["outcomes"][0]["action"], "skipped")

    def test_a_second_run_of_pandda_gets_its_own_receipts(self):
        self.fanout()
        report = self.fanout("--run-job", "99999999-8888-7777-6666-555555555555")
        self.assertEqual(report["outcomes"][0]["action"], "created")
        self.assertEqual(models.Job.objects.filter(project=self.p0).count(), 2, "a rerun is a new receipt")

    def test_plan_is_pure(self):
        plan = plan_fanout(self.tree, json.loads(self.manifest.read_text()), None)
        self.assertEqual(plan.count("created"), 2)
        self.assertEqual(models.Job.objects.count(), 0)
