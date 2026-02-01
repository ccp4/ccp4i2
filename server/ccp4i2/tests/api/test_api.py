"""
API tests for CCP4i2 endpoints.

Converted to pytest fixture-based approach for compatibility with pytest-xdist
parallel test execution. Uses the isolated_test_db fixture from conftest.py.
"""
from pathlib import Path
import json
import pytest
from rest_framework.test import APIClient
from django.core.files.uploadedfile import SimpleUploadedFile

from ccp4i2.core import CCP4Container
from ccp4i2.db.import_i2xml import import_i2xml_from_file, import_ccp4_project_zip
from ccp4i2.db import models

# API URL prefix - all API endpoints are under /api/ccp4i2/
API_PREFIX = "/api/ccp4i2"

# Path to test data - these tests require pre-built project zips
TEST_DATA_DIR = Path(__file__).parent.parent.parent.parent.parent.parent / "test101" / "ProjectZips"
SKIP_REASON = f"Test data not found: {TEST_DATA_DIR}"


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestCCP4i2API:
    """Tests for CCP4i2 API endpoints"""

    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions, test_project_path):
        """Set up test fixtures for each test."""
        self.client = APIClient()
        self.test_project_path = test_project_path
        test_project_path.mkdir(exist_ok=True)
        import_ccp4_project_zip(
            TEST_DATA_DIR / "refmac_gamma_test_0.ccp4_project.zip",
            relocate_path=test_project_path,
        )
        import_i2xml_from_file(
            Path(__file__).parent.parent / "db" / "DATABASE.db.xml",
            relocate_path=test_project_path,
        )

    def test_import_test_dbxml(self):
        assert len(list(models.Project.objects.all())) == 2

    def test_projects(self):
        response = self.client.get(
            f"{API_PREFIX}/projects/", {"username": "john", "password": "smith"}
        )
        project_list = response.json()
        assert project_list[1]["name"] == "MDM2CCP4X"

    def test_project_files(self):
        response = self.client.get(
            f"{API_PREFIX}/projects/2/files/",
        )
        assert len(response.json()) == 24

    def test_project_tags(self):
        response = self.client.get(
            f"{API_PREFIX}/projects/2/tags/",
        )
        assert len(response.json()) == 1

    def test_clone(self):
        response = self.client.post(
            f"{API_PREFIX}/jobs/1/clone/",
        )
        result = response.json()
        # Check key fields - IDs and numbers may vary based on test execution order
        assert "number" in result
        assert int(result["number"]) >= 2  # Cloned job should have number >= 2
        assert result["status"] == 1  # PENDING
        assert result["task_name"] == "prosmart_refmac"

    def test_set_simple_parameter(self):
        # Clone the job first since we can't modify completed jobs
        clone_response = self.client.post(f"{API_PREFIX}/jobs/1/clone/")
        clone = clone_response.json()

        response = self.client.post(
            f"{API_PREFIX}/jobs/{clone['id']}/set_parameter/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(
                {
                    "object_path": "prosmart_refmac.controlParameters.NCYCLES",
                    "value": 20,
                }
            ),
        )
        result = response.json()
        assert result == {
            "success": True,
            "data": {"updated_item": "<NCYCLES>20</NCYCLES>"},
        }

    def test_set_file(self):
        # Clone the job first since we can't modify completed jobs
        clone_response = self.client.post(f"{API_PREFIX}/jobs/1/clone/")
        clone = clone_response.json()

        response = self.client.post(
            f"{API_PREFIX}/jobs/{clone['id']}/set_parameter/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(
                {
                    "object_path": "prosmart_refmac.inputData.XYZIN",
                    "value": {"dbFileId": "AFILEID", "subType": 1},
                }
            ),
        )
        result = response.json()
        assert result == {
            "success": True,
            "data": {"updated_item": "<XYZIN><dbFileId>AFILEID</dbFileId><subType>1</subType></XYZIN>"},
        }

    def test_set_file_null(self):
        # Clone the job first since we can't modify completed jobs
        clone_response = self.client.post(f"{API_PREFIX}/jobs/1/clone/")
        clone = clone_response.json()

        response = self.client.post(
            f"{API_PREFIX}/jobs/{clone['id']}/set_parameter/",
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
            f"{API_PREFIX}/projects/1/preview_file/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(
                {
                    "path": "CCP4_IMPORTED_FILES/gamma_model_1.pdb",
                    "viewer": "coot",
                }
            ),
        )
        result = response.json()
        assert result == {
            "success": True,
            "data": {},
        }

    def test_file_upload(self):
        clone_response = self.client.post(
            f"{API_PREFIX}/jobs/1/clone/",
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
            f"{API_PREFIX}/jobs/{clone['id']}/upload_file_param/", data, format="multipart"
        )
        assert response.json()["data"]["updated_item"]["_value"]["baseName"]["_value"] == "testfile.pdb"

    @pytest.mark.skip(reason="Uploading to empty list index [0] not yet supported - requires list auto-expansion")
    def test_pdb_item_file_upload(self):
        project = models.Project.objects.last()
        create_response = self.client.post(
            f"{API_PREFIX}/projects/{project.id}/create_task/",
            data=json.dumps({"task_name": "phaser_simple"}),
            content_type="application/json; charset=utf-8",
        )
        print(create_response.json())
        # Create a sample file
        file_content = b"Hello, this is a test file."
        test_file = SimpleUploadedFile("testfile.pdb", file_content)

        # Create the data to be sent in the request
        # NOTE: This path references ENSEMBLES[0].pdbItemList[0] which doesn't exist
        # in a freshly created phaser_simple task - the lists are empty
        data = {
            "file": test_file,
            "objectPath": "phaser_simple.inputData.ENSEMBLES[0].pdbItemList[0].structure",
        }
        response = self.client.post(
            f"{API_PREFIX}/jobs/{create_response.json()['data']['new_job']['id']}/upload_file_param/",
            data,
            format="multipart",
        )
        assert response.json()["data"]["updated_item"]["_value"]["baseName"]["_value"] == "testfile.pdb"

    def test_project_upload(self):
        # Create a sample file
        test_file = None
        from django.conf import settings
        try:
            test_project_path = (
                settings.BASE_DIR.parent.parent.parent
                / "test101"
                / "ProjectZips"
                / "refmac_gamma_test_0.ccp4_project.zip"
            )
            test_file = open(test_project_path, "rb")
        except Exception as err:
            print(f"Error opening test file: {err}")
            file_content = b"Hello, this is a test file."
            test_file = SimpleUploadedFile("testfile.pdb", file_content)

        # Create the data to be sent in the request
        data = {"files": [test_file]}
        response = self.client.post(
            f"{API_PREFIX}/projects/import_project/", data, format="multipart"
        )
        print(response)

    def test_mmcif_file_upload(self):
        project = models.Project.objects.last()
        create_response = self.client.post(
            f"{API_PREFIX}/projects/{project.id}/create_task/",
            data=json.dumps({"task_name": "import_merged"}),
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
                f"{API_PREFIX}/jobs/{import_merged_task['new_job']['id']}/upload_file_param/",
                data,
                format="multipart",
            )
            digest_url = f"{API_PREFIX}/jobs/{import_merged_task['new_job']['id']}/digest/?object_path=import_merged.inputData.HKLIN/"
            digest_response = self.client.get(
                digest_url, content_type="application/json; charset=utf-8"
            )
            digest = digest_response.json()["data"]
            assert digest["digest"]["cell"] == {
                "a": 71.45,
                "b": 71.45,
                "c": 104.204,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 120.0,
            }

    def test_digest_file(self):
        file = models.File.objects.first()
        digest_url = f"{API_PREFIX}/files/{file.id}/digest/"
        digest_response = self.client.get(
            digest_url, content_type="application/json; charset=utf-8"
        )
        assert digest_response.json()["data"] == {
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
        }

    def test_upload_to_ProvideAsuContent(self):
        project = models.Project.objects.last()
        create_response = self.client.post(
            f"{API_PREFIX}/projects/{project.id}/create_task/",
            data=json.dumps({"task_name": "ProvideAsuContents"}),
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
                f"{API_PREFIX}/jobs/{ProvideAsuContentsTask['new_job']['id']}/upload_file_param/",
                data,
                format="multipart",
            )
            uploaded_file_uuid = response.json()["data"]["updated_item"]["_value"]["dbFileId"][
                "_value"
            ]
            uploaded_file = models.File.objects.get(uuid=uploaded_file_uuid)
            assert uploaded_file.name == "gamma.pir"
            assert uploaded_file.path.exists()
            digest_url = f"{API_PREFIX}/files/{uploaded_file.id}/digest/"
            digest_response = self.client.get(
                digest_url, content_type="application/json; charset=utf-8"
            )
            print(digest_response.json())
