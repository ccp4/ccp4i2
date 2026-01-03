from pathlib import Path
from shutil import rmtree
import json
import logging
from django.test import Client
from django.conf import settings
from django.test import TestCase, override_settings
from django.core.files.uploadedfile import SimpleUploadedFile
from ccp4i2.core import CCP4Container
from ...db.import_i2xml import import_i2xml_from_file
from ...db.import_i2xml import import_ccp4_project_zip
from ...db import models

logger = logging.getLogger(f"ccp4i2::{__name__}")


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
        import_i2xml_from_file(
            Path(__file__).parent.parent / "db" / "DATABASE.db.xml",
            relocate_path=settings.CCP4I2_PROJECTS_DIR,
        )
        self.client = Client()
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def test_import_test_dbxml(self):
        self.assertEqual(len(list(models.Project.objects.all())), 2)

    def test_projects(self):
        response = self.client.get(
            "/projects/", {"username": "john", "password": "smith"}
        )
        project_list = response.json()
        self.assertEqual(project_list[1]["name"], "MDM2CCP4X")

    def test_project_files(self):
        response = self.client.get(
            "/projects/2/files/",
        )
        self.assertEqual(len(response.json()), 24)

    def test_project_tags(self):
        response = self.client.get(
            "/projects/2/tags/",
        )
        self.assertEqual(len(response.json()), 1)

    def test_clone(self):
        response = self.client.post(
            "/jobs/1/clone/",
        )
        self.assertDictContainsSubset(
            {
                "id": 13,
                "number": "2",
                "title": "Refinement - Refmacat/Refmac5",
                "status": 1,
                "evaluation": 0,
                "comments": "",
                "finish_time": "1970-01-01T00:00:00Z",
                "task_name": "prosmart_refmac",
                "process_id": None,
                "project": 1,
                "parent": None,
            },
            response.json(),
        )

    def test_set_simple_parameter(self):
        response = self.client.post(
            "/jobs/1/set_parameter/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(
                {
                    "object_path": "prosmart_refmac.controlParameters.NCYCLES",
                    "value": 20,
                }
            ),
        )
        result = response.json()
        self.assertDictEqual(
            result,
            {
                "success": True,
                "data": {"updated_item": "<NCYCLES>20</NCYCLES>"},
            },
        )

    def test_set_file(self):
        response = self.client.post(
            "/jobs/1/set_parameter/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(
                {
                    "object_path": "prosmart_refmac.inputData.XYZIN",
                    "value": {"dbFileId": "AFILEID", "subType": 1},
                }
            ),
        )
        result = response.json()
        self.assertDictEqual(
            result,
            {
                "success": True,
                "data": {"updated_item": "<XYZIN><dbFileId>AFILEID</dbFileId><subType>1</subType></XYZIN>"},
            },
        )

    def test_set_file_null(self):
        response = self.client.post(
            "/jobs/1/set_parameter/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(
                {
                    "object_path": "prosmart_refmac.inputData.XYZIN",
                    "value": None,
                }
            ),
        )
        result = response.json()
        print(result)

    def test_preview_file(self):
        response = self.client.post(
            "/projects/1/preview_file/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(
                {
                    "path": "CCP4_IMPORTED_FILES/gamma_model_1.pdb",
                    "viewer": "coot",
                }
            ),
        )
        result = response.json()
        self.assertDictEqual(
            result,
            {
                "success": True,
                "data": {},
            },
        )

    def test_file_upload(self):
        clone_response = self.client.post(
            "/jobs/1/clone/",
        )
        clone = clone_response.json()
        print(clone)
        # Create a sample file
        file_content = b"Hello, this is a test file."
        test_file = SimpleUploadedFile("testfile.pdb", file_content)

        # Create the data to be sent in the request
        data = {
            "file": test_file,
            "objectPath": "prosmart_refmac.inputData.XYZIN",
        }
        response = self.client.post(
            f"/jobs/{clone['id']}/upload_file_param/", data, format="multipart"
        )
        self.assertEqual(
            response.json()["data"]["updated_item"]["_value"]["baseName"]["_value"],
            "testfile.pdb",
        )

    def test_pdb_item_file_upload(self):
        project = models.Project.objects.last()
        create_response = self.client.post(
            f"/projects/{project.id}/create_task/",
            {"task_name": "phaser_simple"},
            content_type="application/json; charset=utf-8",
        )
        print(create_response.json())
        # Create a sample file
        file_content = b"Hello, this is a test file."
        test_file = SimpleUploadedFile("testfile.pdb", file_content)

        # Create the data to be sent in the request
        data = {
            "file": test_file,
            "objectPath": "phaser_simple.inputData.ENSEMBLES[0].pdbItemList[0].structure",
        }
        response = self.client.post(
            f"/jobs/{create_response.json()['data']['new_job']['id']}/upload_file_param/",
            data,
            format="multipart",
        )
        self.assertEqual(
            response.json()["data"]["updated_item"]["_value"]["baseName"]["_value"],
            "testfile.pdb",
        )

    def test_project_upload(self):
        # Create a sample file
        test_file = None
        try:
            test_project_path = (
                settings.BASE_DIR.parent.parent.parent
                / "test101"
                / "ProjectZips"
                / "refmac_gamma_test_0.ccp4_project.zip"
            )
            test_file = open(test_project_path, "rb")
        except Exception as err:
            logger.error(f"Error opening test file: {err}")
            file_content = b"Hello, this is a test file."
            test_file = SimpleUploadedFile("testfile.pdb", file_content)

        # Create the data to be sent in the request
        data = {"files": [test_file]}
        response = self.client.post(
            "/projects/import_project/", data, format="multipart"
        )
        print(response)

    def test_mmcif_file_upload(self):
        project = models.Project.objects.last()
        create_response = self.client.post(
            f"/projects/{project.id}/create_task/",
            {"task_name": "import_merged"},
            content_type="application/json; charset=utf-8",
        )
        import_merged_task = create_response.json()["data"]
        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "4hg7-sf.cif"
        )
        print(mmcif_path, mmcif_path.exists())
        with open(mmcif_path, "rb") as mmcif_file:
            # Create the data to be sent in the request
            data = {
                "file": mmcif_file,
                "objectPath": "import_merged.inputData.HKLIN",
            }
            response = self.client.post(
                f"/jobs/{import_merged_task['new_job']['id']}/upload_file_param/",
                data,
                format="multipart",
            )
            digest_url = f"/jobs/{import_merged_task['new_job']['id']}/digest/?object_path=import_merged.inputData.HKLIN/"
            digest_response = self.client.get(
                digest_url, content_type="application/json; charset=utf-8"
            )
            digest = digest_response.json()["data"]
            self.assertDictEqual(
                digest["digest"]["cell"],
                {
                    "a": 71.45,
                    "b": 71.45,
                    "c": 104.204,
                    "alpha": 90.0,
                    "beta": 90.0,
                    "gamma": 120.0,
                },
            )

    def test_digest_file(self):
        file = models.File.objects.first()
        digest_url = f"/files/{file.id}/digest/"
        digest_response = self.client.get(
            digest_url, content_type="application/json; charset=utf-8"
        )
        self.assertDictEqual(
            digest_response.json()["data"],
            {
                "sequences": {
                    "A": "MIPSITAYSKNGLKIEFTFERSNTNPSVTVITIQASNSTELDMTDFVFQAAVPKTFQLQLLSPSSSVVPAFNTGTITQVIKVLNPQKQQLRMRIKLTYNHKGSAMQDLAEVNNFPPQSWQ"
                },
                "composition": {
                    "chains": ["A"],
                    "peptides": ["A"],
                    "nucleics": [],
                    "solventChains": [],
                    "monomers": [],
                    "nresSolvent": 0,
                    "moleculeType": ["PROTEIN"],
                    "containsHydrogen": False,
                },
            },
        )

    def test_upload_to_ProvideAsuContent(self):
        project = models.Project.objects.last()
        create_response = self.client.post(
            f"/projects/{project.id}/create_task/",
            {"task_name": "ProvideAsuContents"},
            content_type="application/json; charset=utf-8",
        )
        ProvideAsuContentsTask = create_response.json()["data"]
        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "gamma"
            / "gamma.pir"
        )
        print(mmcif_path, mmcif_path.exists())
        with open(mmcif_path, "rb") as mmcif_file:
            # Create the data to be sent in the request
            data = {
                "file": mmcif_file,
                "objectPath": "ProvideAsuContents.inputData.ASU_CONTENT[0].source",
            }
            response = self.client.post(
                f"/jobs/{ProvideAsuContentsTask['new_job']['id']}/upload_file_param/",
                data,
                format="multipart",
            )
            uploaded_file_uuid = response.json()["data"]["updated_item"]["_value"]["dbFileId"][
                "_value"
            ]
            uploaded_file = models.File.objects.get(uuid=uploaded_file_uuid)
            self.assertEqual(uploaded_file.name, "gamma.pir")
            self.assertTrue(uploaded_file.path.exists())
            digest_url = f"/files/{uploaded_file.id}/digest/"
            digest_response = self.client.get(
                digest_url, content_type="application/json; charset=utf-8"
            )
            print(digest_response.json())
