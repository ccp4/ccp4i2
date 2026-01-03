from pathlib import Path
from shutil import rmtree
from django.test import TestCase, override_settings
from django.conf import settings
from xml.etree import ElementTree as ET
from ...db.models import Job
from ...db.import_i2xml import import_ccp4_project_zip
from ...lib.utils.reporting.i2_report import (
    get_report_job_info,
    generate_job_report,
)
from ...lib.utils.jobs.get_container import get_job_container
from ...db.ccp4i2_django_projects_manager import CCP4i2DjangoProjectsManager


from ...db.ccp4i2_django_dbapi import CCP4i2DjangoDbApi


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent.parent.parent.parent
    / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class CCP4i2TestCase(TestCase):
    def setUp(self):
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir()
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "refmac_gamma_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def test_status_labels(self):
        the_job = Job.objects.get(id=1)
        self.assertEqual("Finished", Job.Status(the_job.status).label)

    def test_get_job_container(self):
        the_job = Job.objects.get(id=1)
        container = get_job_container(the_job)
        self.assertListEqual(
            container.dataOrder(),
            [
                "inputData",
                "outputData",
                "controlParameters",
                "prosmartProtein",
                "prosmartNucleicAcid",
                "libg",
                "platonyzer",
                "guiAdmin",
            ],
        )

    def test_get_report_job_info(self):
        the_job = Job.objects.get(id=1)
        result = get_report_job_info(the_job.uuid)
        self.assertEqual(result["finishtime"], 1594563149.18)

    def test_generate_job_report(self):
        the_job = Job.objects.get(id=1)
        self.assertTrue(isinstance(generate_job_report(the_job), ET.Element))

    def test_ccp4_db(self):
        a = CCP4i2DjangoDbApi()
        self.assertTrue(isinstance(a, CCP4i2DjangoDbApi))

    def test_decorator(self):
        def test_fn():
            return 1

        test_fn()
        self.assertEqual(1, 1)

    def test_ccp4_projects_manager(self):
        _ = CCP4i2DjangoProjectsManager()
        self.assertEqual(1, 1)
