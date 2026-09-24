"""Fan-in (design note section 3.1): a pending pandda_campaign job in a
campaign's parent project is filled from the members that have a finished
dimple job, each dataset carrying its own project uuid; members without one
are named with a reason; a second fill adds nothing; a project that is not
a campaign parent says so. Exercised the way the interface reaches it: as
plugin methods through the generic object_method helper."""
import xml.etree.ElementTree as ET
from pathlib import Path
from shutil import rmtree

from django.test import TestCase, override_settings

from ...db import models
from ...lib.utils.helpers.object_method import object_method
from ...lib.utils.jobs.create import create_job
from ...lib.utils.jobs.pandda_fanin import campaign_candidates, fill_datasets_from_campaign

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_FANIN_TEST_DIR"
PARENT = models.ProjectGroupMembership.MembershipType.PARENT
MEMBER = models.ProjectGroupMembership.MembershipType.MEMBER


def make_project(name):
    directory = PROJECTS_DIR / name
    (directory / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    return models.Project.objects.create(name=name, directory=str(directory))


def register(job, param_name, name, mime):
    file_type, _ = models.FileType.objects.get_or_create(name=mime, defaults={"description": mime})
    (job.directory / name).write_text("x\n")
    return models.File.objects.create(name=name, directory=models.File.Directory.JOB_DIR,
                                      type=file_type, job=job, job_param_name=param_name)


def finished_dimple(project, number="1", with_files=True, registered=True):
    job = models.Job.objects.create(project=project, number=number, task_name="i2Dimple",
                                    title="dimple", status=models.Job.Status.FINISHED)
    job.directory.mkdir(parents=True, exist_ok=True)
    if with_files and registered:
        register(job, "XYZOUT", "final.pdb", "chemical/x-pdb")
        register(job, "COMPLETE_MTZ", "final.mtz", "application/CCP4-mtz")
    elif with_files:
        (job.directory / "final.pdb").write_text("END\n")
        (job.directory / "final.mtz").write_bytes(b"MTZ ")
    return job


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class FanInTest(TestCase):
    def setUp(self):
        PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)
        self.parent = make_project("BAZ2B-ref")
        self.m1, self.m2, self.m3 = (make_project(n) for n in ("x425", "x427", "x428"))
        self.group = models.ProjectGroup.objects.create(name="BAZ2B", type="fragment_set")
        for project, kind in ((self.parent, PARENT), (self.m1, MEMBER), (self.m2, MEMBER), (self.m3, MEMBER)):
            models.ProjectGroupMembership.objects.create(group=self.group, project=project, type=kind)
        self.d1 = finished_dimple(self.m1)
        acedrg = models.Job.objects.create(project=self.m1, number="2", task_name="LidiaAcedrgNew",
                                           title="acedrg", status=models.Job.Status.FINISHED)
        acedrg.directory.mkdir(parents=True, exist_ok=True)
        # named after the ligand, as acedrg does, and registered as a list item
        self.dict_row = register(acedrg, "DICTOUT_LIST[0]", "MZ0.cif", "application/refmac-dictionary")
        self.d2 = finished_dimple(self.m2, with_files=False)   # ran, but outputs gone
        # m3 has no dimple at all
        job_uuid = create_job(projectId=str(self.parent.uuid), taskName="pandda_campaign", title="PanDDA")
        self.job = models.Job.objects.get(uuid=job_uuid)

    def datasets(self):
        root = ET.parse(Path(self.job.directory) / "input_params.xml").getroot()
        return root.findall("ccp4i2_body/inputData/DATASETS/CPanddaDataset")

    def test_candidates_say_what_and_why(self):
        preview = campaign_candidates(self.job)
        self.assertIsNone(preview["reason"])
        self.assertEqual([c["name"] for c in preview["campaigns"]], ["BAZ2B"])
        self.assertEqual([c["label"] for c in preview["candidates"]], ["x425"])
        c = preview["candidates"][0]
        self.assertEqual(c["project_uuid"], str(self.m1.uuid))
        self.assertEqual(c["source_job_uuid"], str(self.d1.uuid))
        self.assertTrue(c["dict"].endswith("MZ0.cif"))
        self.assertEqual(set(c["refs"]), {"xyzin", "hklin", "dict"})
        self.assertEqual(c["refs"]["dict"]["dbFileId"], self.dict_row.uuid.hex)
        self.assertEqual(c["refs"]["xyzin"]["project"], self.m1.uuid.hex, "the member's project, not the job's")
        self.assertFalse(c["listed"])
        reasons = {s["project"]: s["reason"] for s in preview["skipped"]}
        self.assertIn("no finished dimple", reasons["x428"])
        self.assertIn("final.pdb/final.mtz", reasons["x427"])

    def test_fill_writes_datasets_with_their_own_project(self):
        result = fill_datasets_from_campaign(self.job)
        self.assertTrue(result.success, result.error)
        self.assertEqual(result.data["added"], ["x425"])
        items = self.datasets()
        self.assertEqual(len(items), 1)
        item = items[0]
        self.assertEqual(item.findtext("DTAG"), "x425")
        self.assertEqual(item.findtext("PROJECT_UUID"), str(self.m1.uuid))
        self.assertEqual(item.findtext("SOURCE_JOB_UUID"), str(self.d1.uuid))
        # Registered outputs are referenced by identity: what the picker
        # shows, what the run records a use of, what the manifest keys on.
        self.assertEqual(item.findtext("XYZIN/baseName"), "final.pdb")
        self.assertEqual(item.findtext("XYZIN/relPath"), "CCP4_JOBS/job_1")
        self.assertEqual(item.findtext("XYZIN/project").replace("-", ""), self.m1.uuid.hex)
        self.assertTrue(item.findtext("XYZIN/dbFileId"))
        self.assertEqual(item.findtext("HKLIN/baseName"), "final.mtz")
        self.assertEqual(item.findtext("DICT/baseName"), "MZ0.cif")
        self.assertEqual(item.findtext("DICT/dbFileId").replace("-", ""), self.dict_row.uuid.hex)

    def test_filling_twice_adds_nothing(self):
        fill_datasets_from_campaign(self.job)
        result = fill_datasets_from_campaign(self.job)
        self.assertTrue(result.success)
        self.assertEqual(result.data["added"], [])
        self.assertTrue(result.data["candidates"][0]["listed"])
        self.assertEqual(len(self.datasets()), 1)

    def test_unregistered_dimple_outputs_fall_back_to_the_path(self):
        m4 = make_project("x429")
        models.ProjectGroupMembership.objects.create(group=self.group, project=m4, type=MEMBER)
        finished_dimple(m4, registered=False)
        preview = campaign_candidates(self.job)
        c = next(c for c in preview["candidates"] if c["label"] == "x429")
        self.assertEqual(c["refs"], {})
        fill_datasets_from_campaign(self.job)
        item = next(i for i in self.datasets() if i.findtext("DTAG") == "x429")
        self.assertTrue(item.findtext("XYZIN/baseName").endswith("x429/CCP4_JOBS/job_1/final.pdb"))

    def test_a_member_project_is_not_a_parent(self):
        job_uuid = create_job(projectId=str(self.m1.uuid), taskName="pandda_campaign", title="PanDDA")
        job = models.Job.objects.get(uuid=job_uuid)
        preview = campaign_candidates(job)
        self.assertIn("not the parent", preview["reason"])
        self.assertFalse(fill_datasets_from_campaign(job).success)

    def test_a_finished_job_cannot_be_filled(self):
        self.job.status = models.Job.Status.FINISHED
        self.job.save()
        self.assertIn("only a pending job", campaign_candidates(self.job)["reason"])

    def test_the_plugin_methods_reach_the_same_code(self):
        """What the interface calls: object_method on the plugin itself."""
        preview = object_method(self.job, "pandda_campaign", "campaignCandidates")
        self.assertEqual([c["label"] for c in preview["candidates"]], ["x425"])
        filled = object_method(self.job, "pandda_campaign", "fillDatasetsFromCampaign")
        self.assertTrue(filled["success"], filled)
        self.assertEqual(filled["data"]["added"], ["x425"])
        self.assertEqual(len(self.datasets()), 1)
