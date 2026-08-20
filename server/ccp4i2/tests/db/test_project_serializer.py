"""How ProjectSerializer treats an update as opposed to a create.

Everything a project's edit page can change goes through this serializer, and
what it does on a create is not what it should do on an update: a create is
about to have a directory made for it under CCP4I2_PROJECTS_DIR and must not
take a name already in use, while an update already has a directory somewhere
and its own name is not a name already in use as far as it is concerned.
"""

from pathlib import Path
from shutil import rmtree

from django.conf import settings
from django.test import TestCase, override_settings

from ...api.serializers import ProjectSerializer
from ...db.models import Project


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class ProjectSerializerTestCase(TestCase):
    def setUp(self):
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir(exist_ok=True)
        self.project = self.create_project("MDM2")
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def create_project(self, name: str) -> Project:
        # "__default__" is what the client sends for "put it in the usual place".
        serializer = ProjectSerializer(data={"name": name, "directory": "__default__"})
        self.assertTrue(serializer.is_valid(), serializer.errors)
        return serializer.save()

    def test_update_may_send_back_its_own_name(self):
        """An edit form sends every field, changed or not."""
        serializer = ProjectSerializer(
            instance=self.project,
            data={"name": "MDM2", "description": "Ligand soaks, 2026 campaign"},
            partial=True,
        )
        self.assertTrue(serializer.is_valid(), serializer.errors)
        project = serializer.save()
        self.assertEqual(project.name, "MDM2")
        self.assertEqual(project.description, "Ligand soaks, 2026 campaign")

    def test_update_still_refuses_a_name_another_project_holds(self):
        self.create_project("MDM4")
        serializer = ProjectSerializer(
            instance=self.project, data={"name": "MDM4"}, partial=True
        )
        self.assertFalse(serializer.is_valid())
        self.assertIn("name", serializer.errors)

    def test_rename_leaves_the_project_where_it_is(self):
        """The directory is named for the project it was created as.

        Renaming does not move it -- relocating a project is what the move
        endpoint is for, and it has files to rewrite as well as bytes to move.
        """
        directory = self.project.directory
        serializer = ProjectSerializer(
            instance=self.project, data={"name": "MDM2_2026"}, partial=True
        )
        self.assertTrue(serializer.is_valid(), serializer.errors)
        project = serializer.save()
        self.assertEqual(project.name, "MDM2_2026")
        self.assertEqual(project.directory, directory)

    def test_create_still_refuses_a_name_already_taken(self):
        serializer = ProjectSerializer(data={"name": "MDM2", "directory": "__default__"})
        self.assertFalse(serializer.is_valid())
        self.assertIn("name", serializer.errors)
