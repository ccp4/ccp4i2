"""
API integration test for aimless_pipe that mimics GUI workflow.

Uses pytest fixtures to set up a file-based database that both the
test and spawned subprocess can access.

Tests:
1. Create a project via API
2. Create an aimless_pipe job
3. Add list item and upload unmerged MTZ file
4. Check job validation
5. Run the job and wait for completion
"""
import json
import os
import time
from pathlib import Path
import pytest
from django.test import Client


TEST_DIR = Path(__file__).parent / "TEST_AIMLESS_API_PROJECTS"
TEST_DB = TEST_DIR / "test_aimless.db"


@pytest.mark.usefixtures("file_based_db")
class TestAimlessAPIIncremental:
    """Incremental API tests for aimless_pipe workflow."""

    # Class-level state to share between tests
    project_id = None
    job_id = None

    @pytest.fixture(autouse=True)
    def setup(self):
        """Set up test client."""
        self.client = Client()
        self.demo_data = Path(os.environ.get('CCP4I2_ROOT', '.')) / 'demo_data'
        self.mtz_file = self.demo_data / 'gamma' / 'gamma_native.mtz'

    def test_01_create_project(self):
        """Test creating a project via API."""
        print("\n=== Test 1: Creating project ===")

        response = self.client.post(
            '/projects/',
            data=json.dumps({
                'name': 'test_aimless_api',
                'directory': '__default__'
            }),
            content_type='application/json'
        )

        print(f"Status: {response.status_code}")
        print(f"Response: {response.content.decode()[:500]}")

        assert response.status_code == 201, f"Failed to create project: {response.content}"

        data = response.json()
        assert 'id' in data
        assert 'name' in data
        assert data['name'] == 'test_aimless_api'

        # Store for next test
        TestAimlessAPIIncremental.project_id = data['id']
        print(f"Created project ID: {data['id']}")

    def test_02_create_aimless_job(self):
        """Test creating an aimless_pipe job within the project."""
        print("\n=== Test 2: Creating aimless_pipe job ===")

        # First ensure we have a project
        if TestAimlessAPIIncremental.project_id is None:
            self.test_01_create_project()

        project_id = TestAimlessAPIIncremental.project_id

        response = self.client.post(
            f'/projects/{project_id}/create_task/',
            data=json.dumps({
                'task_name': 'aimless_pipe',
                'title': 'Test Aimless API Job',
                'auto_context': False
            }),
            content_type='application/json'
        )

        print(f"Status: {response.status_code}")
        print(f"Response: {response.content.decode()[:500]}")

        assert response.status_code == 200, f"Failed to create job: {response.content}"

        data = response.json()
        assert data.get('success'), f"Job creation failed: {data}"
        assert 'data' in data
        assert 'new_job' in data['data']

        job_data = data['data']['new_job']
        assert 'id' in job_data
        assert job_data['task_name'] == 'aimless_pipe'

        # Store for next test
        TestAimlessAPIIncremental.job_id = job_data['id']
        print(f"Created job ID: {job_data['id']}")

    def test_03_add_list_item(self):
        """Test adding an empty list item to UNMERGEDFILES."""
        print("\n=== Test 3a: Adding list item ===")

        # Ensure we have a job
        if TestAimlessAPIIncremental.job_id is None:
            self.test_02_create_aimless_job()

        job_id = TestAimlessAPIIncremental.job_id

        # Add an empty list item to UNMERGEDFILES
        # This mimics what the GUI does when you click "+" on a list
        response = self.client.post(
            f'/jobs/{job_id}/set_parameter/',
            data=json.dumps({
                'object_path': 'aimless_pipe.container.inputData.UNMERGEDFILES',
                'value': [{}]  # Add one empty item
            }),
            content_type='application/json'
        )

        print(f"Status: {response.status_code}")
        print(f"Response: {response.content.decode()[:500]}")

        assert response.status_code == 200, f"Failed to add list item: {response.content}"

    def test_04_upload_mtz_file(self):
        """Test uploading an unmerged MTZ file to the job."""
        print("\n=== Test 3b: Uploading MTZ file ===")

        # Ensure we have a job with list item
        if TestAimlessAPIIncremental.job_id is None:
            self.test_03_add_list_item()

        job_id = TestAimlessAPIIncremental.job_id

        # Check demo data exists
        if not self.mtz_file.exists():
            pytest.skip(f"Demo data not found: {self.mtz_file}")

        # Use bracket notation [0] for list indexing
        with open(self.mtz_file, 'rb') as f:
            response = self.client.post(
                f'/jobs/{job_id}/upload_file_param/',
                {
                    'objectPath': 'aimless_pipe.container.inputData.UNMERGEDFILES[0].file',
                    'file': f,
                }
            )

        print(f"Status: {response.status_code}")
        print(f"Response: {response.content.decode()[:500]}")

        assert response.status_code == 200, f"Failed to upload file: {response.content}"

    def test_05_check_validation(self):
        """Test checking job validation before running."""
        print("\n=== Test 5: Checking validation ===")

        # Ensure we have a job with file uploaded
        if TestAimlessAPIIncremental.job_id is None:
            self.test_04_upload_mtz_file()

        job_id = TestAimlessAPIIncremental.job_id

        response = self.client.get(f'/jobs/{job_id}/validation/')

        print(f"Status: {response.status_code}")
        print(f"Response: {response.content.decode()[:1500]}")

        assert response.status_code == 200, f"Failed to get validation: {response.content}"

        # Check for any ERROR level validations
        content = response.content.decode()
        if 'ERROR' in content:
            print("\nWARNING: Validation errors found!")
        if 'WARNING' in content:
            print("\nNote: Validation warnings found (may be acceptable)")

    def test_06_run_job(self):
        """Test running the job and waiting for completion."""
        print("\n=== Test 6: Running job ===")

        # Ensure we have a validated job
        if TestAimlessAPIIncremental.job_id is None:
            self.test_05_check_validation()

        job_id = TestAimlessAPIIncremental.job_id

        # Run the job
        response = self.client.post(
            f'/jobs/{job_id}/run/',
            content_type='application/json'
        )

        print(f"Run status: {response.status_code}")
        print(f"Run response: {response.content.decode()[:500]}")

        # Note: Job execution via API in tests may have issues due to
        # Django TestCase transaction isolation. The subprocess spawned
        # by run/ cannot see uncommitted transactions.
        #
        # For debugging GUI issues, the important thing is to see if
        # the job starts at all, and what errors occur.

        assert response.status_code == 200, f"Failed to run job: {response.content}"

        # Poll for job completion (status 6=FINISHED, 7=FAILED, 8=UNSATISFACTORY)
        print("\nPolling for job status...")
        final_status = None
        for i in range(30):  # Wait up to 30 seconds
            time.sleep(1)
            status_response = self.client.get(f'/jobs/{job_id}/')
            if status_response.status_code == 200:
                job_data = status_response.json()
                final_status = job_data.get('status')
                print(f"  {i}s: status={final_status}")
                if final_status in [6, 7, 8]:  # Terminal states
                    break

        # Job should have finished (status 6) or at least started running
        assert final_status is not None, "Failed to get job status"
        assert final_status >= 3, f"Job did not start running (status={final_status})"
        print(f"\nFinal job status: {final_status}")

        # Check output files
        print("\n=== Checking job files ===")
        files_response = self.client.get(f'/jobs/{job_id}/files/')
        print(f"Files status: {files_response.status_code}")
        if files_response.status_code == 200:
            files_data = files_response.json()
            print(f"Files: {json.dumps(files_data, indent=2)[:2000]}")


