"""
API tests for context-aware parameter setting endpoints.

Tests the set_parameter endpoint to verify that:
1. Parameters can be set through HTTP POST
2. Database synchronization occurs automatically via CPluginScript context
3. Various parameter types are handled correctly (int, str, bool, files)
4. Error cases are handled appropriately

Converted to pytest fixture-based approach for compatibility with pytest-xdist
parallel test execution. Uses the isolated_test_db fixture from conftest.py.
"""

import json
import uuid
import pytest
from pathlib import Path
from rest_framework.test import APIClient

from ccp4i2.lib.utils.parameters.get_param import get_parameter
from ccp4i2.db.import_i2xml import import_ccp4_project_zip
from ccp4i2.db import models

# API URL prefix - all API endpoints are under /api/ccp4i2/
API_PREFIX = "/api/ccp4i2"

# Path to test data - these tests require pre-built project zips
TEST_DATA_DIR = Path(__file__).parent.parent.parent.parent.parent.parent / "test101" / "ProjectZips"
SKIP_REASON = f"Test data not found: {TEST_DATA_DIR}"


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestParameterSettingAPI:
    """Tests for parameter setting API endpoints with context-aware database sync"""

    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions, test_project_path):
        """Set up test project and client"""
        self.client = APIClient()
        self.test_project_path = test_project_path
        test_project_path.mkdir(parents=True, exist_ok=True)

        import_ccp4_project_zip(
            TEST_DATA_DIR / "refmac_gamma_test_0.ccp4_project.zip",
            relocate_path=test_project_path,
        )

        # Get the project to create a new job in
        project = models.Project.objects.first()
        assert project is not None, "Should have at least one project"

        # IMPORTANT: Create a fresh PENDING job with proper CData wrappers
        # This matches i2run workflow: load plugin from .def.xml, set params, run
        job_number = str(models.Job.objects.filter(project=project).count() + 1)

        # Create job in database (directory is computed from project + number)
        self.job = models.Job.objects.create(
            uuid=uuid.uuid4(),
            project=project,
            number=job_number,
            task_name="prosmart_refmac",
            title="Test Parameter Setting",
            status=models.Job.Status.PENDING
        )

        # Create job directory
        job_dir = self.job.directory
        job_dir.mkdir(parents=True, exist_ok=True)

        # Verify job is PENDING and ready for parameter setting
        assert self.job.status == models.Job.Status.PENDING

    def test_set_parameter_simple_integer(self):
        """Test setting a simple integer parameter via API"""
        # Test data: set NCYCLES parameter (in controlParameters, not inputData)
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES",
            "value": 5
        }

        # POST to set_parameter endpoint
        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        # Verify response
        assert response.status_code == 200
        response_data = response.json()
        print(f"\n📤 API Response: {response_data}")
        assert response_data.get("success")
        assert "data" in response_data

        # DEBUG: Check if input_params.xml contains NCYCLES
        input_params_file = self.job.directory / "input_params.xml"
        print(f"\n🔍 DEBUG: Checking for NCYCLES in {input_params_file}")
        if input_params_file.exists():
            xml_content = input_params_file.read_text()
            print("✅ input_params.xml exists, searching for NCYCLES...")
            if "NCYCLES" in xml_content:
                print("✅ NCYCLES found in input_params.xml")
                # Extract the value to see what was written
                import re
                match = re.search(r'NCYCLES[^>]*>(\d+)', xml_content)
                if match:
                    print(f"✅ NCYCLES value in XML: {match.group(1)}")
            else:
                print("❌ NCYCLES NOT found in input_params.xml")
                print(f"XML content (first 2000 chars):\n{xml_content[:2000]}")
        else:
            print(f"❌ input_params.xml does NOT exist at {input_params_file}")

        # Verify parameter was actually set
        result = get_parameter(self.job, "prosmart_refmac.controlParameters.NCYCLES")
        assert result.success, f"get_parameter failed: {result.error if not result.success else 'no error'}"
        assert result.data["value"] == 5

    def test_set_parameter_boolean(self):
        """Test setting a boolean parameter via API"""
        # Test data: set ADD_WATERS parameter
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.ADD_WATERS",
            "value": False
        }

        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        assert response.status_code == 200
        response_data = response.json()
        assert response_data.get("success")

        # Verify parameter was set
        result = get_parameter(self.job, "prosmart_refmac.controlParameters.ADD_WATERS")
        assert result.success
        assert result.data["value"] == False

    def test_set_parameter_string(self):
        """Test setting a string parameter via API"""
        # Test data: set job title
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.TITLE",
            "value": "API Test Job"
        }

        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        assert response.status_code == 200
        response_data = response.json()
        assert response_data.get("success")

        # Verify parameter was set
        result = get_parameter(self.job, "prosmart_refmac.controlParameters.TITLE")
        assert result.success
        assert result.data["value"] == "API Test Job"

    def test_set_parameter_invalid_path(self):
        """Test error handling for invalid parameter path"""
        request_data = {
            "object_path": "prosmart_refmac.inputData.NONEXISTENT_PARAM",
            "value": 42
        }

        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        # Should return error status
        assert response.status_code == 400
        response_data = response.json()
        assert not response_data.get("success")
        assert "error" in response_data

    def test_set_parameter_invalid_job_id(self):
        """Test error handling for non-existent job"""
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES",
            "value": 5
        }

        # Use non-existent job ID (999999)
        response = self.client.post(
            f"{API_PREFIX}/jobs/999999/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        # Should return 404
        assert response.status_code == 404
        response_data = response.json()
        assert not response_data.get("success")
        assert "Job not found" in response_data.get("error", "")

    def test_set_parameter_missing_required_fields(self):
        """Test error handling for missing required fields"""
        # Missing 'value' field
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES"
        }

        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        # Should return error (400 or 500 depending on validation)
        assert response.status_code in [400, 500]
        response_data = response.json()
        assert not response_data.get("success")

    def test_set_parameter_database_sync(self):
        """Test that parameter setting triggers database synchronization"""
        # Set a parameter
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES",
            "value": 10
        }

        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        assert response.status_code == 200

        # Verify job was updated in database by checking input_params.xml exists
        input_params_file = self.job.directory / "input_params.xml"
        assert input_params_file.exists(), \
            f"input_params.xml should exist after parameter setting at {input_params_file}"

        # Read and verify the XML contains our parameter
        xml_content = input_params_file.read_text()
        # The XML should contain some representation of NCYCLES
        assert "NCYCLES" in xml_content, "NCYCLES parameter should be in saved XML"

    def test_set_parameter_multiple_sequential(self):
        """Test setting multiple parameters sequentially"""
        parameters = [
            ("prosmart_refmac.controlParameters.NCYCLES", 5),
            ("prosmart_refmac.controlParameters.ADD_WATERS", True),
            ("prosmart_refmac.controlParameters.VALIDATE_MOLPROBITY", False),
        ]

        for param_path, param_value in parameters:
            request_data = {
                "object_path": param_path,
                "value": param_value
            }

            response = self.client.post(
                f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
                data=json.dumps(request_data),
                content_type="application/json"
            )

            assert response.status_code == 200, f"Failed to set parameter {param_path}"

            # Verify each parameter was set
            result = get_parameter(self.job, param_path)
            assert result.success, f"Failed to retrieve parameter {param_path}"
            assert result.data["value"] == param_value, f"Parameter {param_path} has wrong value"

    def test_set_parameter_null_value(self):
        """Test setting a parameter to null/None - should fail for typed fields like CInt"""
        # First set a parameter
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES",
            "value": 5
        }
        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )
        assert response.status_code == 200

        # Setting None on a CInt should fail - CInt doesn't accept None values
        # To unset a parameter, use the unSet() method on the CData object
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES",
            "value": None
        }
        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        # CInt cannot be set to None - should return error
        assert response.status_code == 400
        response_data = response.json()
        assert not response_data.get("success")

    def test_set_parameter_type_coercion(self):
        """Test that string values are coerced to correct types"""
        # Set integer parameter with string value
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES",
            "value": "7"  # String instead of int
        }

        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        assert response.status_code == 200

        # Verify it was coerced to integer
        result = get_parameter(self.job, "prosmart_refmac.controlParameters.NCYCLES")
        assert result.success
        # Should be integer, not string
        assert isinstance(result.data["value"], int)
        assert result.data["value"] == 7

    def test_set_parameter_with_skip_first(self):
        """Test parameter path handling with skip_first behavior"""
        # Some paths may include task name as first element
        # The API should handle this correctly
        request_data = {
            "object_path": "prosmart_refmac.controlParameters.NCYCLES",
            "value": 3
        }

        response = self.client.post(
            f"{API_PREFIX}/jobs/{self.job.id}/set_parameter/",
            data=json.dumps(request_data),
            content_type="application/json"
        )

        assert response.status_code == 200
        response_data = response.json()
        assert response_data.get("success")
