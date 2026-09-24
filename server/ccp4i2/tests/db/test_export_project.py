from pathlib import Path
from shutil import rmtree
from xml.etree import ElementTree as ET
from glob import glob
from django.test import TestCase, override_settings

from .external_data import (
    PROJECTS_SCRATCH_DIR,
    TEST_ZIPS_DIR,
    requires_project_zips,
)
from django.conf import settings
from ...db.models import Project
from ...db.import_i2xml import import_i2xml_from_file, import_ccp4_project_zip
from ...db.export_project import generate_project_xml_tree, export_project_to_xml


@requires_project_zips
@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_SCRATCH_DIR)
class CCP4i2TestCase(TestCase):
    def setUp(self):
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir(parents=True, exist_ok=True)
        import_ccp4_project_zip(
            TEST_ZIPS_DIR
            / "refmac_gamma_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        import_ccp4_project_zip(
            TEST_ZIPS_DIR
            / "aimless_gamma_native_test_1.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        import_ccp4_project_zip(
            TEST_ZIPS_DIR
            / "parrot_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR, ignore_errors=True)
        return super().tearDown()

    def test_export_project(self):
        project: Project = Project.objects.get(name="refmac_gamma_test_0")
        a = generate_project_xml_tree(project)
        ET.indent(a, space="  ")
        print(ET.tostring(a, encoding="unicode"))
