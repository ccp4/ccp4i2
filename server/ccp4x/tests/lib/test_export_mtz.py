from pathlib import Path
from shutil import rmtree
from xml.etree import ElementTree as ET
from django.test import TestCase, override_settings
from django.conf import settings
from ...db.models import Job, File, FileImport
from ...db.import_i2xml import import_ccp4_project_zip

from ...lib.utils.plugins.get_plugin import get_job_plugin
from ...lib.utils.formats.mtz import mtz_as_dict
from ...lib.utils.parameters.unset_output_data import unset_output_data
from ...lib.utils.containers.remove_defaults import (
    remove_container_default_values,
)
from ...lib.utils.containers.find_objects import find_objects
from ...lib.utils.parameters.load_xml import load_nested_xml
from ...lib.utils.containers.validate import validate_container
from ...lib.utils.jobs.clone import clone_job
from ...lib.utils.jobs.create import create_job
from ...lib.utils.containers.json_for_container import json_for_job_container
from ...lib.utils.navigation.task_tree import get_task_tree
from ...lib.utils.reporting.i2_report import get_report_job_info
from ...lib.utils.formats.gemmi_split_mtz import gemmi_split_mtz
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
