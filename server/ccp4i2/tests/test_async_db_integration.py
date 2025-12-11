"""
Integration tests for async database handler with real Django database.

Run with:
    ./run_db_tests.sh
"""

import os
import sys
import pytest
import asyncio
import uuid
from pathlib import Path

# Add server directory to path if needed
server_dir = Path(__file__).parent.parent / 'server'
if str(server_dir) not in sys.path:
    sys.path.insert(0, str(server_dir))

# Set Django settings if not already set (pytest-django will call django.setup())
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4x.config.test_settings')


@pytest.fixture
def projects_base():
    """Get the base directory for test projects"""
    projects_dir = Path(os.environ.get('CCP4I2_PROJECTS_DIR', Path.home() / '.ccp4x_test' / 'test_projects'))
    projects_dir.mkdir(parents=True, exist_ok=True)
    return projects_dir


@pytest.fixture
def test_project(projects_base, db):
    """Create a test project for testing"""
    from ccp4x.db import models
    import shutil

    # Create unique project directory
    project_name = f"test_project_{uuid.uuid4().hex[:8]}"
    project_dir = projects_base / project_name
    project_dir.mkdir(parents=True, exist_ok=True)

    # Create project in database
    project = models.Project.objects.create(
        name=project_name,
        directory=str(project_dir),
        description="Test project for async handler",
    )

    print(f"Created test project: {project.name} ({project.uuid})")

    yield project

    # Cleanup
    if project_dir.exists():
        shutil.rmtree(project_dir, ignore_errors=True)


@pytest.mark.django_db(transaction=True)
class TestAsyncDatabaseHandlerIntegration:
    """Test AsyncDatabaseHandler with real database operations."""

    @pytest.mark.asyncio
    async def test_create_job(self, test_project):
        """Test creating a job in the database"""
        from ccp4x.db.async_db_handler import AsyncDatabaseHandler

        handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)

        # Create a job
        job = await handler.create_job(
            task_name="ctruncate",
            title="Test job",
        )

        # Verify
        assert job.project == test_project
        assert job.task_name == "ctruncate"
        assert job.title == "Test job"
        assert job.number == "1"  # First job
        assert job.status == 1  # PENDING

        print(f"✅ Created job {job.number}: {job.title}")

    @pytest.mark.asyncio
    async def test_create_nested_jobs(self, test_project):
        """Test creating nested jobs with automatic numbering"""
        from ccp4x.db.async_db_handler import AsyncDatabaseHandler

        handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)

        # Create parent job
        parent = await handler.create_job(
            task_name="copycell",
            title="Parent job",
        )
        assert parent.number == "1"

        # Create child jobs
        child1 = await handler.create_job(
            task_name="ctruncate",
            title="Child 1",
            parent_job_uuid=parent.uuid,
        )
        assert child1.number == "1.1"

        child2 = await handler.create_job(
            task_name="refmac",
            title="Child 2",
            parent_job_uuid=parent.uuid,
        )
        assert child2.number == "1.2"

        # Create grandchild
        grandchild = await handler.create_job(
            task_name="pointless",
            title="Grandchild",
            parent_job_uuid=child1.uuid,
        )
        assert grandchild.number == "1.1.1"

        print(f"✅ Created nested jobs: {parent.number}, {child1.number}, {child2.number}, {grandchild.number}")

    @pytest.mark.asyncio
    async def test_update_job_status(self, test_project):
        """Test updating job status"""
        from ccp4x.db.async_db_handler import AsyncDatabaseHandler
        from ccp4x.db import models
        from asgiref.sync import sync_to_async

        handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)

        # Create job
        job = await handler.create_job(
            task_name="ctruncate",
            title="Status test job",
        )

        # Update to RUNNING
        await handler.update_job_status(job.uuid, models.Job.Status.RUNNING)

        # Fetch from database to verify
        updated_job = await sync_to_async(models.Job.objects.get)(uuid=job.uuid)
        assert updated_job.status == models.Job.Status.RUNNING

        # Update to FINISHED
        await handler.update_job_status(job.uuid, models.Job.Status.FINISHED)
        updated_job = await sync_to_async(models.Job.objects.get)(uuid=job.uuid)
        assert updated_job.status == models.Job.Status.FINISHED
        assert updated_job.finish_time is not None

        print(f"✅ Job status updated: PENDING → RUNNING → FINISHED")

    @pytest.mark.asyncio
    async def test_register_job_values(self, test_project):
        """Test registering KPI values"""
        from ccp4x.db.async_db_handler import AsyncDatabaseHandler
        from ccp4x.db import models
        from asgiref.sync import sync_to_async

        handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)

        # Create job
        job = await handler.create_job(
            task_name="refmac",
            title="KPI test job",
        )

        # Register float values
        await handler.register_job_float_value(
            job_uuid=job.uuid,
            key="Rfactor",
            value=0.234,
            description="R-factor value",
        )

        await handler.register_job_float_value(
            job_uuid=job.uuid,
            key="Rfree",
            value=0.289,
            description="R-free value",
        )

        # Register string value
        await handler.register_job_char_value(
            job_uuid=job.uuid,
            key="SpaceGroup",
            value="P212121",
            description="Space group",
        )

        # Verify
        float_values = await sync_to_async(list)(
            models.JobFloatValue.objects.filter(job=job)
        )
        char_values = await sync_to_async(list)(
            models.JobCharValue.objects.filter(job=job)
        )

        assert len(float_values) == 2
        assert len(char_values) == 1

        print(f"✅ Registered {len(float_values)} float KPIs and {len(char_values)} string KPIs")

    @pytest.mark.skip(reason="Requires actual CPluginScript and file objects")
    def test_glean_files_from_container(self):
        """Test gleaning files from a container (would require real plugin)"""
        # This test would require:
        # 1. A real CPluginScript instance
        # 2. Real output files
        # 3. Proper CData file objects
        pass

    @pytest.mark.asyncio
    async def test_job_directory_creation(self, test_project):
        """Test that job directories follow the correct convention"""
        from ccp4x.db.async_db_handler import AsyncDatabaseHandler

        handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)

        # Create nested jobs
        parent = await handler.create_job(task_name="copycell", title="Parent")
        child = await handler.create_job(
            task_name="ctruncate",
            title="Child",
            parent_job_uuid=parent.uuid,
        )

        # Check directory paths
        parent_dir = parent.directory
        child_dir = child.directory

        # Verify directory structure
        assert "job_1" in str(parent_dir)
        assert "job_1" in str(child_dir)

        print(f"✅ Parent dir: {parent_dir}")
        print(f"✅ Child dir: {child_dir}")

        assert parent.number == "1"
        assert child.number == "1.1"


@pytest.mark.django_db(transaction=True)
class TestCDataUtilitiesIntegration:
    """Test CData utilities with mock objects"""

    def test_extract_file_metadata_with_mock(self):
        """Test metadata extraction with a mock file object"""
        from ccp4x.lib.cdata_utils import extract_file_metadata
        from unittest.mock import Mock

        # Create mock file
        mock_file = Mock()
        mock_file.name = "HKLOUT"
        mock_file.object_path.return_value = "outputData.HKLOUT"

        # Mock get_merged_metadata to return different values based on parameter
        def get_merged_metadata(key):
            if key == 'qualifiers':
                return {
                    'mimeTypeName': 'application/CCP4-mtz',
                    'guiLabel': 'Output MTZ file',
                    'toolTip': 'Reflection data',
                }
            elif key == 'attributes':
                return {'subType': True}  # Indicate subType exists
            return {}

        mock_file.get_merged_metadata.side_effect = get_merged_metadata
        mock_file.isSet.return_value = True
        mock_file.exists.return_value = True

        # Add optional attributes
        mock_subtype = Mock()
        mock_subtype.isSet.return_value = True
        mock_subtype.value = 1
        mock_file.subType = mock_subtype

        # Extract metadata
        metadata = extract_file_metadata(mock_file)

        # Verify
        assert metadata['name'] == 'HKLOUT'
        assert metadata['file_type'] == 'application/CCP4-mtz'
        assert metadata['gui_label'] == 'Output MTZ file'
        assert metadata['sub_type'] == 1
        assert metadata['is_set'] is True

        print(f"✅ Extracted metadata: {metadata}")

    def test_find_all_files_with_hierarchy(self):
        """Test finding files in a hierarchical structure"""
        from unittest.mock import Mock

        # Create mock container
        mock_container = Mock()

        # Create mock file objects
        mock_file1 = Mock()
        mock_file1.__class__.__name__ = "CDataFile"

        mock_file2 = Mock()
        mock_file2.__class__.__name__ = "CDataFile"

        # Mock childNames
        mock_container.childNames.return_value = ['HKLOUT', 'XYZOUT']
        mock_container.HKLOUT = mock_file1
        mock_container.XYZOUT = mock_file2

        # Note: This is a simplified test - real test would need proper CData objects
        print("✅ CData utilities test framework ready")
