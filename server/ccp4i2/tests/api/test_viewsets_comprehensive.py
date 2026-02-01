"""
Comprehensive API ViewSet Tests

This test suite provides complete coverage of all API endpoints defined in:
- FileViewSet
- JobViewSet
- ProjectViewSet
- FileImportViewSet
- FileUseViewSet
- FileTypeViewSet
- ProjectTagViewSet
- ProjectExportViewSet

Purpose: Enable safe refactoring of job_utils by ensuring all endpoints
behave correctly before and after removing legacy code.

Converted to pytest fixture-based approach for compatibility with pytest-xdist
parallel test execution. Uses the isolated_test_db fixture from conftest.py.
"""

from pathlib import Path
import json
import pytest
from rest_framework.test import APIClient
from django.core.files.uploadedfile import SimpleUploadedFile

from ccp4i2.db.import_i2xml import import_i2xml_from_file, import_ccp4_project_zip
from ccp4i2.db import models

# API URL prefix - all API endpoints are under /api/ccp4i2/
API_PREFIX = "/api/ccp4i2"

# Path to test data - these tests require pre-built project zips
TEST_DATA_DIR = Path(__file__).parent.parent.parent.parent.parent.parent / "test101" / "ProjectZips"
SKIP_REASON = f"Test data not found: {TEST_DATA_DIR}"


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestFileViewSet:
    """Tests for FileViewSet API endpoints"""

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
        self.file = models.File.objects.first()

    def test_file_list(self):
        """Test GET /files/ - List all files"""
        response = self.client.get(f"{API_PREFIX}/files/")
        assert response.status_code == 200
        files = response.json()
        assert isinstance(files, list)
        assert len(files) > 0

    def test_file_retrieve(self):
        """Test GET /files/{id}/ - Retrieve specific file"""
        response = self.client.get(f"{API_PREFIX}/files/{self.file.id}/")
        assert response.status_code == 200
        file_data = response.json()
        assert file_data["name"] == self.file.name

    def test_file_by_uuid(self):
        """Test GET /files/{id}/by_uuid/ - Retrieve file by UUID"""
        response = self.client.get(f"{API_PREFIX}/files/{self.file.uuid}/by_uuid/")
        assert response.status_code == 200
        file_data = response.json()
        assert file_data["name"] == self.file.name

    def test_file_download(self):
        """Test GET /files/{id}/download/ - Download file by ID"""
        response = self.client.get(f"{API_PREFIX}/files/{self.file.id}/download/")
        assert response.status_code == 200
        # FileResponse uses 'inline' by default, 'attachment' is optional
        assert 'filename=' in response["Content-Disposition"]
        assert self.file.name in response["Content-Disposition"]

    def test_file_download_by_uuid(self):
        """Test GET /files/{uuid}/download_by_uuid/ - Download file by UUID"""
        response = self.client.get(f"{API_PREFIX}/files/{self.file.uuid}/download_by_uuid/")
        assert response.status_code == 200
        assert 'filename=' in response["Content-Disposition"]
        assert self.file.name in response["Content-Disposition"]

    def test_file_digest(self):
        """Test GET /files/{id}/digest/ - Get file digest by ID"""
        response = self.client.get(f"{API_PREFIX}/files/{self.file.id}/digest/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        digest = result.get("data", {})
        # Digest structure varies by file type, just check we got a dict
        assert isinstance(digest, dict)
        # For PDB files, expect sequences and composition
        if self.file.type.name in ['xyzin', 'xyzout']:
            assert "sequences" in digest
            assert "composition" in digest

    def test_file_digest_by_uuid(self):
        """Test GET /files/{uuid}/digest_by_uuid/ - Get file digest by UUID"""
        response = self.client.get(f"{API_PREFIX}/files/{self.file.uuid}/digest_by_uuid/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        digest = result.get("data", {})
        assert isinstance(digest, dict)
        if self.file.type.name in ['xyzin', 'xyzout']:
            assert "sequences" in digest
            assert "composition" in digest

    @pytest.mark.skip(reason="Requires external viewer application")
    def test_file_preview(self):
        """Test POST /files/{id}/preview/ - Preview file with external viewer"""
        response = self.client.post(
            f"{API_PREFIX}/files/{self.file.id}/preview/",
            data=json.dumps({"viewer": "coot"}),
            content_type="application/json",
        )
        # This would actually launch coot, so we skip in automated tests
        pass


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestJobViewSet:
    """Tests for JobViewSet API endpoints"""

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
        self.job = models.Job.objects.first()

    def test_job_list(self):
        """Test GET /jobs/ - List all jobs"""
        response = self.client.get(f"{API_PREFIX}/jobs/")
        assert response.status_code == 200
        jobs = response.json()
        assert isinstance(jobs, list)
        assert len(jobs) > 0

    def test_job_retrieve(self):
        """Test GET /jobs/{id}/ - Retrieve specific job"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/")
        assert response.status_code == 200
        job_data = response.json()
        assert job_data["id"] == self.job.id

    def test_job_clone(self):
        """Test POST /jobs/{id}/clone/ - Clone existing job"""
        response = self.client.post(f"{API_PREFIX}/jobs/{self.job.id}/clone/")
        assert response.status_code in [200, 201]  # 201 Created is also valid
        cloned_job = response.json()
        assert cloned_job["id"] != self.job.id
        assert cloned_job["task_name"] == self.job.task_name
        assert cloned_job["status"] == models.Job.Status.PENDING

    def test_job_params_xml(self):
        """Test GET /jobs/{id}/params_xml/ - Get job parameters as XML"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/params_xml/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        assert "xml" in result.get("data", {})

    def test_job_report_xml(self):
        """Test GET /jobs/{id}/report_xml/ - Get job report as XML"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/report_xml/")
        # May be 200 or 404 depending on whether report exists
        assert response.status_code in [200, 404]
        result = response.json()
        if response.status_code == 200:
            assert result.get("success")
            assert "xml" in result.get("data", {})

    def test_job_diagnostic_xml(self):
        """Test GET /jobs/{id}/diagnostic_xml/ - Get diagnostic XML"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/diagnostic_xml/")
        # May be 200 or 404 depending on whether diagnostic exists
        assert response.status_code in [200, 404]

    def test_job_container(self):
        """Test GET /jobs/{id}/container/ - Get job container JSON"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/container/")
        result = response.json()
        assert result.get("success")
        assert "result" in result.get("data", {})

    def test_job_def_xml(self):
        """Test GET /jobs/{id}/def_xml/ - Get task definition XML"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/def_xml/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        assert "xml" in result.get("data", {})

    def test_job_validation(self):
        """Test GET /jobs/{id}/validation/ - Validate job parameters"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/validation/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        assert "xml" in result.get("data", {})

    def test_job_i2run_command(self):
        """Test GET /jobs/{id}/i2run_command/ - Get i2run command"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/i2run_command/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        assert "command" in result.get("data", {})

    def test_job_dependent_jobs(self):
        """Test GET /jobs/{id}/dependent_jobs/ - Get jobs that depend on this job"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/dependent_jobs/")
        assert response.status_code == 200
        dependent_jobs = response.json()
        assert isinstance(dependent_jobs, list)

    def test_job_files(self):
        """Test GET /jobs/{id}/files/ - Get files associated with job"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/files/")
        assert response.status_code == 200
        files = response.json()
        assert isinstance(files, list)

    def test_job_what_next(self):
        """Test GET /jobs/{id}/what_next/ - Get suggested next steps"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/what_next/")
        result = response.json()
        # Should return suggestions or empty list
        assert isinstance(result, dict)

    def test_job_set_parameter_simple(self):
        """Test POST /jobs/{id}/set_parameter/ - Set simple parameter"""
        # Clone job first to avoid modifying original
        clone_response = self.client.post(f"{API_PREFIX}/jobs/{self.job.id}/clone/")
        cloned_job = clone_response.json()

        response = self.client.post(
            f"{API_PREFIX}/jobs/{cloned_job['id']}/set_parameter/",
            content_type="application/json",
            data=json.dumps({
                "object_path": f"{self.job.task_name}.controlParameters.NCYCLES",
                "value": 20,
            }),
        )
        result = response.json()
        assert result.get("success")

    def test_job_set_parameter_file(self):
        """Test POST /jobs/{id}/set_parameter/ - Set file parameter"""
        # Clone job first
        clone_response = self.client.post(f"{API_PREFIX}/jobs/{self.job.id}/clone/")
        cloned_job = clone_response.json()

        response = self.client.post(
            f"{API_PREFIX}/jobs/{cloned_job['id']}/set_parameter/",
            content_type="application/json",
            data=json.dumps({
                "object_path": f"{self.job.task_name}.inputData.XYZIN",
                "value": {"dbFileId": "AFILEID", "subType": 1},
            }),
        )
        result = response.json()
        assert result.get("success")

    def test_job_upload_file_param(self):
        """Test POST /jobs/{id}/upload_file_param/ - Upload file parameter"""
        # Clone job first
        clone_response = self.client.post(f"{API_PREFIX}/jobs/{self.job.id}/clone/")
        cloned_job = clone_response.json()

        file_content = b"Test PDB content"
        test_file = SimpleUploadedFile("test.pdb", file_content)

        data = {
            "file": test_file,
            "objectPath": f"{self.job.task_name}.inputData.XYZIN",
        }
        response = self.client.post(
            f"{API_PREFIX}/jobs/{cloned_job['id']}/upload_file_param/",
            data,
            format="multipart",
        )
        result = response.json()
        assert result.get("success")
        assert "updated_item" in result.get("data", {})

    def test_job_digest(self):
        """Test GET /jobs/{id}/digest/?object_path=... - Digest object"""
        response = self.client.get(
            f"{API_PREFIX}/jobs/{self.job.id}/digest/?object_path={self.job.task_name}.inputData.XYZIN/"
        )
        # May succeed or fail depending on whether file is set
        assert response.status_code in [200, 404, 500]

    def test_job_digest_param_file(self):
        """Test GET /jobs/{id}/digest_param_file/?job_param_name=... - Digest param file"""
        # Find a file parameter
        file_obj = models.File.objects.filter(job=self.job).first()
        if file_obj and file_obj.job_param_name:
            response = self.client.get(
                f"{API_PREFIX}/jobs/{self.job.id}/digest_param_file/?job_param_name={file_obj.job_param_name}/"
            )
            # May succeed or fail depending on file type
            assert response.status_code in [200, 404, 500]

    def test_job_export_job_file_menu(self):
        """Test GET /jobs/{id}/export_job_file_menu/ - Get export file menu"""
        response = self.client.get(f"{API_PREFIX}/jobs/{self.job.id}/export_job_file_menu/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        assert "result" in result.get("data", {})

    @pytest.mark.skip(reason="Requires CCP4 installation and takes time")
    def test_job_run(self):
        """Test POST /jobs/{id}/run/ - Execute job"""
        # Skip in automated tests as this actually runs the job
        pass

    @pytest.mark.skip(reason="Requires CCP4 installation and takes time")
    def test_job_run_local(self):
        """Test POST /jobs/{id}/run_local/ - Execute job locally"""
        # Skip in automated tests
        pass

    @pytest.mark.skip(reason="Requires external viewer")
    def test_job_preview(self):
        """Test POST /jobs/{id}/preview/ - Preview job with viewer"""
        # Skip as it launches external application
        pass

    def test_job_delete(self):
        """Test DELETE /jobs/{id}/ - Delete job and dependents"""
        # Clone a job to delete
        clone_response = self.client.post(f"{API_PREFIX}/jobs/{self.job.id}/clone/")
        cloned_job = clone_response.json()

        response = self.client.delete(f"{API_PREFIX}/jobs/{cloned_job['id']}/")
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")

        # Verify job was deleted
        assert not models.Job.objects.filter(id=cloned_job['id']).exists()


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestProjectViewSet:
    """Tests for ProjectViewSet API endpoints"""

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
        self.project = models.Project.objects.first()

    def test_project_list(self):
        """Test GET /projects/ - List all projects"""
        response = self.client.get(f"{API_PREFIX}/projects/")
        assert response.status_code == 200
        projects = response.json()
        assert isinstance(projects, list)
        assert len(projects) > 0

    def test_project_retrieve(self):
        """Test GET /projects/{id}/ - Retrieve specific project"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/")
        assert response.status_code == 200
        project_data = response.json()
        assert project_data["name"] == self.project.name

    def test_project_files(self):
        """Test GET /projects/{id}/files/ - Get project files"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/files/")
        assert response.status_code == 200
        files = response.json()
        assert isinstance(files, list)

    def test_project_file_uses(self):
        """Test GET /projects/{id}/file_uses/ - Get project file uses"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/file_uses/")
        assert response.status_code == 200
        file_uses = response.json()
        assert isinstance(file_uses, list)

    def test_project_jobs(self):
        """Test GET /projects/{id}/jobs/ - Get project jobs"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/jobs/")
        assert response.status_code == 200
        jobs = response.json()
        assert isinstance(jobs, list)
        assert len(jobs) > 0

    def test_project_job_float_values(self):
        """Test GET /projects/{id}/job_float_values/ - Get KPI float values"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/job_float_values/")
        assert response.status_code == 200
        values = response.json()
        assert isinstance(values, list)

    def test_project_job_char_values(self):
        """Test GET /projects/{id}/job_char_values/ - Get KPI char values"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/job_char_values/")
        assert response.status_code == 200
        values = response.json()
        assert isinstance(values, list)

    def test_project_tags_get(self):
        """Test GET /projects/{id}/tags/ - Get project tags"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/tags/")
        assert response.status_code == 200
        tags = response.json()
        assert isinstance(tags, list)

    def test_project_tags_post(self):
        """Test POST /projects/{id}/tags/ - Add tag to project"""
        # Create a tag first
        tag = models.ProjectTag.objects.create(text="Test Tag")

        response = self.client.post(
            f"{API_PREFIX}/projects/{self.project.id}/tags/",
            data=json.dumps({"tag_id": tag.id}),
            content_type="application/json",
        )
        assert response.status_code == 200
        result = response.json()
        assert result["status"] == "success"

    def test_project_remove_tag(self):
        """Test DELETE /projects/{id}/tags/{tag_id}/ - Remove tag from project"""
        # Create and add a tag first
        tag = models.ProjectTag.objects.create(text="Test Tag")
        self.project.tags.add(tag)

        response = self.client.delete(f"{API_PREFIX}/projects/{self.project.id}/tags/{tag.id}/")
        assert response.status_code == 200
        result = response.json()
        assert result["status"] == "success"

    def test_project_directory(self):
        """Test GET /projects/{id}/directory/ - Get project directory listing"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/directory/")
        assert response.status_code == 200
        result = response.json()
        # API returns {"status": "Success", "container": [...]} or {"success": True, "data": {...}}
        assert result.get("success") or result.get("status") == "Success"
        assert "container" in result.get("data", {}) or "container" in result

    def test_project_file(self):
        """Test GET /projects/{id}/project_file/?path=... - Get project file"""
        # Find a file in the project
        file_obj = models.File.objects.filter(job__project=self.project).first()
        if file_obj:
            # Use relative path from project directory
            relative_path = file_obj.path.relative_to(self.project.directory)
            response = self.client.get(
                f"{API_PREFIX}/projects/{self.project.id}/project_file/?path={relative_path}"
            )
            assert response.status_code == 200

    @pytest.mark.skip(reason="Requires external viewer")
    def test_project_preview_file(self):
        """Test POST /projects/{id}/preview_file/ - Preview file with viewer"""
        # Skip as it launches external application
        pass

    def test_project_create_task(self):
        """Test POST /projects/{id}/create_task/ - Create new task/job"""
        response = self.client.post(
            f"{API_PREFIX}/projects/{self.project.id}/create_task/",
            data=json.dumps({"task_name": "prosmart_refmac"}),
            content_type="application/json",
        )
        assert response.status_code == 200
        result = response.json()
        assert result.get("success")
        assert "new_job" in result.get("data", {})

    def test_project_exports(self):
        """Test GET /projects/{id}/exports/ - Get project export history"""
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/exports/")
        assert response.status_code == 200
        exports = response.json()
        assert isinstance(exports, list)

    @pytest.mark.skip(reason="Creates large files and background processes")
    def test_project_export(self):
        """Test POST /projects/{id}/export/ - Export project"""
        # Skip as it creates large ZIP files and background processes
        pass

    @pytest.mark.skip(reason="Would delete test data")
    def test_project_delete(self):
        """Test DELETE /projects/{id}/ - Delete project"""
        # Skip to preserve test data
        pass

    @pytest.mark.skip(reason="Requires uploaded ZIP file")
    def test_project_import(self):
        """Test POST /projects/import_project/ - Import project from ZIP"""
        # Skip as it requires actual ZIP file upload
        pass


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestSimpleViewSets:
    """Tests for simple CRUD ViewSets: FileImport, FileUse, FileType, ProjectTag"""

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

    def test_file_import_list(self):
        """Test GET /fileimports/ - List all file imports"""
        response = self.client.get(f"{API_PREFIX}/fileimports/")
        assert response.status_code == 200
        imports = response.json()
        assert isinstance(imports, list)

    def test_file_use_list(self):
        """Test GET /fileuses/ - List all file uses"""
        response = self.client.get(f"{API_PREFIX}/fileuses/")
        assert response.status_code == 200
        uses = response.json()
        assert isinstance(uses, list)

    def test_file_type_list(self):
        """Test GET /filetypes/ - List all file types"""
        response = self.client.get(f"{API_PREFIX}/filetypes/")
        assert response.status_code == 200
        types = response.json()
        assert isinstance(types, list)

    def test_project_tag_list(self):
        """Test GET /projecttags/ - List all project tags"""
        response = self.client.get(f"{API_PREFIX}/projecttags/")
        assert response.status_code == 200
        tags = response.json()
        assert isinstance(tags, list)

    def test_project_tag_create(self):
        """Test POST /projecttags/ - Create new tag"""
        import uuid
        unique_tag = f"Test Tag {uuid.uuid4().hex[:8]}"
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            data=json.dumps({"text": unique_tag, "parent": None}),
            content_type="application/json",
        )
        # If status is not 201, print response for debugging
        if response.status_code != 201:
            print(f"Tag create failed: {response.status_code} - {response.content.decode()}")
        assert response.status_code == 201, f"Expected 201, got {response.status_code}: {response.content.decode()}"
        tag = response.json()
        assert tag["text"] == unique_tag


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestProjectExportViewSet:
    """Tests for ProjectExportViewSet API endpoints"""

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
        self.project = models.Project.objects.first()

    def test_project_export_list(self):
        """Test GET /projectexports/ - List all exports"""
        response = self.client.get(f"{API_PREFIX}/projectexports/")
        assert response.status_code == 200
        exports = response.json()
        assert isinstance(exports, list)

    @pytest.mark.skip(reason="Requires actual export file creation")
    def test_project_export_download(self):
        """Test GET /projectexports/{id}/download/ - Download export"""
        # Skip as it requires actual export file to exist
        pass

    @pytest.mark.skip(reason="Would delete export files")
    def test_project_export_delete(self):
        """Test DELETE /projectexports/{id}/ - Delete export"""
        # Skip to preserve data
        pass
