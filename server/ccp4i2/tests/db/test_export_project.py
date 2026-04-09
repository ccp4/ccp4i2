from pathlib import Path
from shutil import rmtree
from xml.etree import ElementTree as ET
from glob import glob
from django.test import TestCase, override_settings
from django.conf import settings
from ...db.models import Project
from ...db.import_i2xml import import_i2xml_from_file, import_ccp4_project_zip
from ...db.export_project import generate_project_xml_tree, export_project_to_xml


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent / "CCP4I2_TEST_PROJECT_DIRECTORY"
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
        rmtree(settings.CCP4I2_PROJECTS_DIR, ignore_errors=True)
        return super().tearDown()

    def test_export_project(self):
        project: Project = Project.objects.get(name="refmac_gamma_test_0")
        a = generate_project_xml_tree(project)
        ET.indent(a, space="  ")
        print(ET.tostring(a, encoding="unicode"))
