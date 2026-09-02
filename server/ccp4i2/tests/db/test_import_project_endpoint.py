"""The upload endpoint has to look inside the zip before accepting it.

``import_project`` dispatches the real work with ``--detach``, so nothing the
child discovers can reach the person who pressed the button: a zip the importer
cannot read dies in a subprocess while the upload reports success and no project
appears. Whatever is worth telling the user therefore has to be found *before*
the dispatch, which is what these tests pin down.
"""

import io
import tempfile
import zipfile
from pathlib import Path
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.test import TestCase, override_settings
from rest_framework.test import APIClient

from ...db.models import Project

PROJECT_XML = """<?xml version='1.0' encoding='ASCII'?>
<ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header>
    <projectName>Demo</projectName>
    <projectId>f48588ee9ca611f1895cacde48001122</projectId>
  </ccp4i2_header>
  <ccp4i2_body>
    <projectTable>
      <project projectid="f48588ee9ca611f1895cacde48001122" projectname="Demo"
               projectdirectory="/elsewhere/CCP4I2_PROJECTS/Demo"/>
    </projectTable>
    <jobTable><job jobid="a" jobnumber="1"/></jobTable>
    <fileTable/>
  </ccp4i2_body>
</ccp4:ccp4i2>
"""

MEDIA_DIR = tempfile.mkdtemp(prefix="ccp4i2-import-endpoint-")


def zip_bytes(members: dict) -> io.BytesIO:
    """An in-memory zip of ``{name: contents}``, ready to upload."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, contents in members.items():
            archive.writestr(name, contents)
    buffer.seek(0)
    buffer.name = "upload.zip"
    return buffer


@override_settings(MEDIA_ROOT=MEDIA_DIR)
class ImportProjectEndpointTest(TestCase):
    def setUp(self):
        self.client = APIClient()
        self.client.force_authenticate(
            get_user_model().objects.create(username="tester")
        )

    def post(self, payload):
        return self.client.post(
            "/api/ccp4i2/projects/import_project/",
            {"files": payload},
            format="multipart",
        )

    @patch("ccp4i2.api.ProjectViewSet.call_command")
    def test_an_exported_zip_is_accepted(self, call_command):
        response = self.post(zip_bytes({"DATABASE.db.xml": PROJECT_XML}))

        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(call_command.call_count, 1)
        described = response.data["data"]["projects"][0]
        self.assertEqual(described["project_name"], "Demo")
        self.assertEqual(described["jobs"], 1)

    @patch("ccp4i2.api.ProjectViewSet.call_command")
    def test_a_hand_rolled_zip_is_accepted_too(self, call_command):
        """``zip -r proj.zip ./PROJ`` wraps the project in a directory."""
        response = self.post(zip_bytes({"PROJ/DATABASE.db.xml": PROJECT_XML}))

        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(call_command.call_count, 1)
        self.assertEqual(response.data["data"]["projects"][0]["project_name"], "Demo")

    @patch("ccp4i2.api.ProjectViewSet.call_command")
    def test_a_zip_of_something_else_is_refused_before_dispatch(self, call_command):
        response = self.post(zip_bytes({"holiday/photo.jpg": "not a project"}))

        self.assertEqual(response.status_code, 400)
        self.assertIn("DATABASE.db.xml", response.data["error"])
        # The point of the pre-flight: nothing was handed to a subprocess that
        # would have failed out of sight.
        call_command.assert_not_called()
        self.assertEqual(Project.objects.count(), 0)

    @patch("ccp4i2.api.ProjectViewSet.call_command")
    def test_a_file_that_is_not_a_zip_is_refused(self, call_command):
        plain = io.BytesIO(b"hello")
        plain.name = "notes.zip"

        response = self.post(plain)

        self.assertEqual(response.status_code, 400)
        call_command.assert_not_called()

    @patch("ccp4i2.api.ProjectViewSet.call_command")
    def test_the_uploaded_file_is_still_saved_for_the_import(self, call_command):
        self.post(zip_bytes({"DATABASE.db.xml": PROJECT_XML}))

        dispatched = Path(call_command.call_args.args[1])
        self.assertTrue(dispatched.is_file())
        self.assertEqual(dispatched.parent, Path(MEDIA_DIR) / "uploaded_files")
