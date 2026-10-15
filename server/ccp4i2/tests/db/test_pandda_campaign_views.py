"""Where a PanDDA run's events are looked at: the campaign the run belongs
to and each dataset's member project with its current-model job, as the
plugin records them in the report's ``<peek>`` so its buttons open the
campaign viewer on the right model. Read through the plugin the way the
report generation reaches it (a plugin that did not run: the manifest comes
from the job directory)."""
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from shutil import rmtree

from django.test import TestCase, override_settings

from ...db import models
from ...lib.utils.helpers.object_method import object_method
from ...lib.utils.jobs.create import create_job
from ...wrappers.pandda_campaign.script.pandda_campaign_report import pandda_campaign_report

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_PANDDA_VIEWS_TEST_DIR"
PARENT = models.ProjectGroupMembership.MembershipType.PARENT
MEMBER = models.ProjectGroupMembership.MembershipType.MEMBER


def make_project(name):
    directory = PROJECTS_DIR / name
    (directory / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    return models.Project.objects.create(name=name, directory=str(directory))


def finished(project, number, task):
    job = models.Job.objects.create(project=project, number=number, task_name=task,
                                    title=task, status=models.Job.Status.FINISHED)
    job.directory.mkdir(parents=True, exist_ok=True)
    return job


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class CampaignViewsTest(TestCase):
    def setUp(self):
        PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)
        self.parent = make_project("BAZ2B-ref")
        self.m1, self.m2 = make_project("x425"), make_project("x427")
        self.group = models.ProjectGroup.objects.create(name="BAZ2B", type="fragment_set")
        for project, kind in ((self.parent, PARENT), (self.m1, MEMBER), (self.m2, MEMBER)):
            models.ProjectGroupMembership.objects.create(group=self.group, project=project, type=kind)
        # m1: dimple, then a refinement; the refinement is the current model.
        finished(self.m1, "1", "i2Dimple")
        self.refined = finished(self.m1, "2", "servalcat_pipe")
        # m2: nothing finished, so no current model.
        job_uuid = create_job(projectId=str(self.parent.uuid), taskName="pandda_campaign", title="PanDDA")
        self.job = models.Job.objects.get(uuid=job_uuid)
        (Path(self.job.directory) / "manifest.json").write_text(json.dumps({
            "manifest_version": 1, "datasets": [
                {"xtal": "xtal-0000", "label": "x425", "project_uuid": str(self.m1.uuid), "files": {}},
                {"xtal": "xtal-0001", "label": "x427", "project_uuid": str(self.m2.uuid), "files": {}},
                {"xtal": "xtal-0002", "label": "gone", "project_uuid": "00000000-0000-0000-0000-000000000000",
                 "files": {}},
            ]}))

    def views(self, job):
        return object_method(job, "pandda_campaign", "_campaign_views", [], {})

    def test_members_are_named_with_their_current_model(self):
        views = self.views(self.job)
        assert views["campaign"] == self.group.id
        assert views["members"] == {
            "xtal-0000": (self.m1.id, self.refined.id),
            "xtal-0001": (self.m2.id, None),
        }

    def test_a_run_outside_a_campaign_names_none(self):
        lone = make_project("lone")
        job_uuid = create_job(projectId=str(lone.uuid), taskName="pandda_campaign", title="PanDDA")
        job = models.Job.objects.get(uuid=job_uuid)
        (Path(job.directory) / "manifest.json").write_text(json.dumps({"manifest_version": 1, "datasets": []}))
        assert self.views(job) == {"campaign": None, "members": {}}

    def test_a_report_drawn_from_an_older_record_looks_the_views_up(self):
        """A run finished before the plugin recorded where events are looked
        at: regenerating its report is enough, the report asks the database."""
        xml = ("<pandda_campaign><state>finished</state><peek total='2' datasets='2'>"
               "<event dtag='xtal-0000' idx='1' score='0.9'/><event dtag='xtal-0001' idx='1'/>"
               "</peek></pandda_campaign>")
        report = pandda_campaign_report(xmlnode=ET.fromstring(xml),
                                        jobInfo={"project_pk": self.parent.id, "jobid": str(self.job.uuid)},
                                        jobStatus="Finished")
        links = [e for e in report.as_data_etree().iter() if e.tag == "CCP4i2ReportMoorhenLink"]
        assert (links[0].get("campaignId"), links[0].get("jobId")) == (str(self.group.id), str(self.refined.id))
        # a member with no current model stays on the project page
        assert links[1].get("campaignId") is None
