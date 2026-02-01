from pathlib import Path
from shutil import rmtree
from django.test import TestCase, override_settings
from django.conf import settings
from ...db.models import Job, File, FileImport
from ...db.import_i2xml import import_ccp4_project_zip

from ...lib.utils.files.export_mtz import export_job_mtz_file, get_source_reflection_file


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
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "aimless_gamma_native_test_1.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "parrot_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )

        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def test_aimless_get_source(self):
        the_job = Job.objects.filter(task_name="prosmart_refmac").first()
        print(File.objects.filter(job=the_job))
        print(FileImport.objects.filter(file__job=the_job))
        get_source_reflection_file(
            jobId=the_job.uuid,
            jobParamNameList=["F_SIGF"],
        )

    def test_export_job_mtz_file(self):
        task_name = "parrot"
        job = Job.objects.filter(task_name=task_name).first()
        if not job:
            raise ValueError(f"No job found with task_name '{task_name}'")
        export_job_mtz_file(job.uuid)
