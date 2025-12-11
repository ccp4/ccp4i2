"""
Tests for set_input_by_context_job functionality.

Tests the context-based input population that occurs when creating a new job
that should inherit inputs from a previous "context" job's outputs.
"""
from pathlib import Path
from shutil import rmtree
import logging

from django.test import TestCase, override_settings
from django.conf import settings

from ...db.models import Job, File, Project
from ...db.import_i2xml import import_i2xml_from_file
from ...lib.utils.parameters.set_input_by_context import set_input_by_context_job
from ...lib.utils.plugins.get_plugin import get_job_plugin
from ...lib.utils.files.get_by_context import get_file_by_job_context
from ...lib.async_create_job import create_job_async
from asgiref.sync import async_to_sync

# Enable logging for debugging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent.parent.parent.parent
    / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class SetInputByContextTestCase(TestCase):
    """Test set_input_by_context_job with imported project data."""

    def setUp(self):
        """Import the test project with aimless_pipe -> prosmart_refmac jobs."""
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir(exist_ok=True)
        import_i2xml_from_file(
            Path(__file__).parent.parent / "db" / "DATABASE.db.xml",
            relocate_path=settings.CCP4I2_PROJECTS_DIR,
        )
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR, ignore_errors=True)
        return super().tearDown()

    def test_project_imported_correctly(self):
        """Verify test data was imported correctly."""
        from ...db.models import FileType

        projects = list(Project.objects.all())
        self.assertEqual(len(projects), 1)
        self.assertEqual(projects[0].name, "MDM2CCP4X")

        jobs = list(Job.objects.filter(parent__isnull=True))
        self.assertEqual(len(jobs), 2)  # aimless_pipe and prosmart_refmac

        # Check FileTypes
        print("\n=== FileTypes ===")
        for ft in FileType.objects.all():
            print(f"  FileType: {ft.name}")

        # Check that aimless_pipe has output files
        aimless_job = Job.objects.get(task_name="aimless_pipe")
        aimless_files = list(File.objects.filter(job=aimless_job))
        print(f"\n=== Aimless job {aimless_job.number} has {len(aimless_files)} files ===")
        for f in aimless_files:
            type_name = f.type.name if f.type else 'None'
            print(f"  - {f.name} (type_id={f.type_id}, type={type_name}, sub_type={f.sub_type}, content={f.content})")

        self.assertGreater(len(aimless_files), 0)

    def test_get_file_by_job_context_finds_mtz(self):
        """Test that get_file_by_job_context can find MTZ files from context job."""
        aimless_job = Job.objects.get(task_name="aimless_pipe")
        project = aimless_job.project

        # The type names in db are MIME types like "application/CCP4-mtz-observed"
        # Try to find observed data MTZ using the actual MIME type
        file_ids = get_file_by_job_context(
            contextJobId=str(aimless_job.uuid),
            fileType="application/CCP4-mtz-observed",  # MIME type, not class name
            subType=1,  # OBSERVED
            contentFlag=1,  # FMEAN
            projectId=str(project.uuid),
        )
        print(f"\nFound {len(file_ids)} application/CCP4-mtz-observed files with subType=1, contentFlag=1")
        self.assertGreater(len(file_ids), 0, "Should find observed data MTZ file")

    def test_get_file_by_job_context_finds_freerflag(self):
        """Test that get_file_by_job_context can find FreeR flag file."""
        aimless_job = Job.objects.get(task_name="aimless_pipe")
        project = aimless_job.project

        # Try to find FreeR flag MTZ using MIME type (database stores MIME types)
        file_ids = get_file_by_job_context(
            contextJobId=str(aimless_job.uuid),
            fileType="application/CCP4-mtz-freerflag",  # MIME type, not class name
            subType=None,
            contentFlag=1,
            projectId=str(project.uuid),
        )
        logger.info(f"Found {len(file_ids)} FreeR flag files")
        self.assertGreater(len(file_ids), 0, "Should find FreeR flag MTZ file")

    def test_get_file_by_job_context_finds_pdb(self):
        """Test that get_file_by_job_context can find PDB coordinate file."""
        aimless_job = Job.objects.get(task_name="aimless_pipe")
        project = aimless_job.project

        # Try to find coordinate file using MIME type (database stores MIME types)
        file_ids = get_file_by_job_context(
            contextJobId=str(aimless_job.uuid),
            fileType="chemical/x-pdb",  # MIME type, not class name
            subType=None,
            contentFlag=None,
            projectId=str(project.uuid),
        )
        logger.info(f"Found {len(file_ids)} PDB files")
        self.assertGreater(len(file_ids), 0, "Should find PDB coordinate file")

    def test_get_job_plugin_loads_plugin(self):
        """Test that get_job_plugin can load a plugin for a job."""
        aimless_job = Job.objects.get(task_name="aimless_pipe")
        plugin = get_job_plugin(aimless_job)

        self.assertIsNotNone(plugin)
        self.assertIsNotNone(plugin.container)
        self.assertIsNotNone(plugin.container.inputData)
        logger.info(f"Plugin container has inputData: {plugin.container.inputData.dataOrder()}")

    def test_set_input_by_context_job_basic(self):
        """Test set_input_by_context_job populates inputs correctly."""
        project = Project.objects.first()
        aimless_job = Job.objects.get(task_name="aimless_pipe")

        # Create a new prosmart_refmac job
        result = async_to_sync(create_job_async)(
            project_uuid=project.uuid,
            task_name="prosmart_refmac",
            title="Test ProSMART Refmac",
            context_job_uuid=aimless_job.uuid,
            auto_context=False,
            save_params=True,
        )

        new_job = Job.objects.get(uuid=result["job_uuid"])
        logger.info(f"Created new job {new_job.number} for task prosmart_refmac")

        # Load the plugin and check if inputs are set
        plugin = get_job_plugin(new_job)
        input_data = plugin.container.inputData

        # Find all input files
        all_files = input_data.find_all_files()
        logger.info(f"New job has {len(all_files)} input file objects")

        set_files = []
        for f in all_files:
            if f.isSet():
                set_files.append(f)
                logger.info(f"  SET: {f.objectName()} = {f.getFullPath()}")
            else:
                logger.info(f"  UNSET: {f.objectName()}")

        # Should have at least some files set from context
        self.assertGreater(
            len(set_files), 0,
            f"Expected some input files to be set from context job. "
            f"Found {len(all_files)} file objects, {len(set_files)} set."
        )

    def test_set_input_by_context_finds_from_previous_job_files(self):
        """Test that files with fromPreviousJob=True qualifier get populated."""
        project = Project.objects.first()
        aimless_job = Job.objects.get(task_name="aimless_pipe")

        # Create a new job without context first
        result = async_to_sync(create_job_async)(
            project_uuid=project.uuid,
            task_name="prosmart_refmac",
            title="Test Without Context",
            context_job_uuid=None,
            auto_context=False,
            save_params=False,
        )

        new_job = Job.objects.get(uuid=result["job_uuid"])
        plugin = result["plugin"]

        # Debug: Check what's in the container
        container = plugin.container
        input_data = container.inputData

        print(f"\n=== Container debug ===")
        print(f"Container type: {type(container)}")
        print(f"Container dataOrder: {container.dataOrder() if hasattr(container, 'dataOrder') else 'N/A'}")
        print(f"InputData type: {type(input_data)}")
        print(f"InputData dataOrder: {input_data.dataOrder() if hasattr(input_data, 'dataOrder') else 'N/A'}")
        print(f"InputData children: {input_data.children() if hasattr(input_data, 'children') else 'N/A'}")

        all_files = input_data.find_all_files()

        print(f"\n=== All {len(all_files)} input files ===")
        for f in all_files:
            class_name = f.__class__.__name__
            obj_name = f.objectName() if hasattr(f, 'objectName') else str(f)
            all_quals = f.qualifiers() if hasattr(f, 'qualifiers') else {}
            from_prev = f.qualifiers("fromPreviousJob") if hasattr(f, 'qualifiers') else None
            print(f"  - {obj_name} ({class_name}): fromPreviousJob={from_prev}, all_quals={all_quals}")

        from_previous_files = [f for f in all_files if f.qualifiers("fromPreviousJob")]

        print(f"\n=== Found {len(from_previous_files)} files with fromPreviousJob=True ===")
        for f in from_previous_files:
            mime_type = f.qualifiers("mimeTypeName")
            sub_type = f.qualifiers("requiredSubType")
            content_flag = f.qualifiers("requiredContentFlag")
            class_name = f.__class__.__name__
            print(f"  - {f.objectName()} ({class_name}): mimeTypeName={mime_type}, subType={sub_type}, contentFlag={content_flag}")

        # Now apply context
        set_input_by_context_job(
            job_id=str(new_job.uuid),
            context_job_id=str(aimless_job.uuid),
            plugin=plugin,
            save_params=False,
        )

        # Check which files got set
        set_count = 0
        for f in from_previous_files:
            if f.isSet():
                set_count += 1
                logger.info(f"  AFTER CONTEXT - SET: {f.objectName()} = {f.getFullPath()}")
            else:
                logger.info(f"  AFTER CONTEXT - UNSET: {f.objectName()}")

        logger.info(f"After context: {set_count}/{len(from_previous_files)} fromPreviousJob files are set")

        # At least some files should be set
        self.assertGreater(
            set_count, 0,
            f"Expected some fromPreviousJob files to be set. Found {len(from_previous_files)} qualifying files."
        )

    def test_set_input_by_context_with_explicit_plugin(self):
        """Test that passing explicit plugin modifies that instance."""
        project = Project.objects.first()
        aimless_job = Job.objects.get(task_name="aimless_pipe")

        # Create job and get plugin
        result = async_to_sync(create_job_async)(
            project_uuid=project.uuid,
            task_name="prosmart_refmac",
            title="Test Explicit Plugin",
            context_job_uuid=None,
            auto_context=False,
            save_params=False,
        )

        new_job = Job.objects.get(uuid=result["job_uuid"])
        plugin = result["plugin"]

        # Count set files before
        input_files = plugin.container.inputData.find_all_files()
        from_previous = [f for f in input_files if f.qualifiers("fromPreviousJob")]
        set_before = sum(1 for f in from_previous if f.isSet())

        # Apply context with explicit plugin
        set_input_by_context_job(
            job_id=str(new_job.uuid),
            context_job_id=str(aimless_job.uuid),
            plugin=plugin,
            save_params=False,
        )

        # Count set files after - should be same plugin instance
        set_after = sum(1 for f in from_previous if f.isSet())

        logger.info(f"Before context: {set_before} set, After context: {set_after} set")
        self.assertGreater(set_after, set_before, "Context should set more files")
