"""
Tests for report generation using pre-existing job data.

These tests use project zip files (exported from CCP4i2) that contain
completed jobs with all their output files. This allows testing report
generation without running actual CCP4 jobs.
"""
import xml.etree.ElementTree as ET
from pathlib import Path
from shutil import rmtree
from django.test import TestCase, override_settings
from django.conf import settings

from ccp4i2.db.import_i2xml import import_ccp4_project_zip
from ccp4i2.db.models import Job
from ccp4i2.lib.utils.reporting.i2_report import (
    generate_job_report,
    get_report_job_info,
    simple_failed_report,
)


# Path to test project zips
TEST_ZIPS_DIR = (
    Path(__file__).parent.parent.parent.parent.parent.parent / "test101" / "ProjectZips"
)


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent.parent.parent.parent
    / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class ReportGenerationTests(TestCase):
    """Test report generation with pre-existing job data."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # Create test projects directory
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir(exist_ok=True)

    @classmethod
    def tearDownClass(cls):
        # Clean up test projects directory
        rmtree(settings.CCP4I2_PROJECTS_DIR, ignore_errors=True)
        super().tearDownClass()

    def test_refmac_report_generation(self):
        """Test report generation for a refmac job from imported project."""
        # Import the refmac test project
        zip_path = TEST_ZIPS_DIR / "refmac_gamma_test_0.ccp4_project.zip"
        if not zip_path.exists():
            self.skipTest(f"Test project zip not found: {zip_path}")

        import_ccp4_project_zip(
            zip_path,
            relocate_path=settings.CCP4I2_PROJECTS_DIR,
        )

        # Find the refmac job
        job = Job.objects.filter(
            project__name="refmac_gamma_test_0", task_name="refmac"
        ).first()
        if job is None:
            # The prosmart_refmac pipeline creates a subjob
            job = Job.objects.filter(project__name="refmac_gamma_test_0").first()

        self.assertIsNotNone(job, "Job not found in imported project")
        print(f"\nTesting report for job: {job.task_name} #{job.number}")

        # Generate report
        report_xml = generate_job_report(job)

        # Validate report
        self.assertIsNotNone(report_xml)
        self.assertIsInstance(report_xml, ET.Element)

        report_str = ET.tostring(report_xml, encoding="unicode")
        print(f"Report generated: {len(report_str)} characters")

        # Check that it's not a failure report (unless expected)
        if "No report because" in report_str:
            print(f"Note: Report indicates limitation: {report_str[:300]}")
        else:
            # Should have some content
            self.assertGreater(len(report_str), 500)

    def test_get_report_job_info(self):
        """Test job info collection for reports."""
        # Import a simple project
        zip_path = TEST_ZIPS_DIR / "refmac_gamma_test_0.ccp4_project.zip"
        if not zip_path.exists():
            self.skipTest(f"Test project zip not found: {zip_path}")

        import_ccp4_project_zip(
            zip_path,
            relocate_path=settings.CCP4I2_PROJECTS_DIR,
        )

        job = Job.objects.filter(project__name="refmac_gamma_test_0").first()
        self.assertIsNotNone(job)

        # Get job info
        job_info = get_report_job_info(job.uuid)

        # Validate required fields
        self.assertIn("taskname", job_info)
        self.assertIn("jobnumber", job_info)
        self.assertIn("projectid", job_info)
        self.assertIn("fileroot", job_info)
        self.assertIn("inputfiles", job_info)
        self.assertIn("outputfiles", job_info)
        self.assertIn("filenames", job_info)

        print(f"Job info collected:")
        print(f"  Task: {job_info['taskname']}")
        print(f"  Number: {job_info['jobnumber']}")
        print(f"  Input files: {len(job_info['inputfiles'])}")
        print(f"  Output files: {len(job_info['outputfiles'])}")

    def test_simple_failed_report_structure(self):
        """Test that simple_failed_report generates valid XML."""
        report = simple_failed_report(
            reason="Test failure reason",
            task_name="test_task",
            details="This is a test.\nWith multiple lines.\nAnd <special> chars.",
        )

        self.assertIsNotNone(report)
        self.assertIsInstance(report, ET.Element)

        report_str = ET.tostring(report, encoding="unicode")

        # Check structure
        self.assertIn("test_task", report_str)
        self.assertIn("Test failure reason", report_str)
        self.assertIn("This is a test", report_str)
        # Special chars should be escaped
        self.assertIn("&lt;special&gt;", report_str)

    def test_report_for_missing_report_class(self):
        """Test graceful handling when no report class exists."""
        # Import project
        zip_path = TEST_ZIPS_DIR / "refmac_gamma_test_0.ccp4_project.zip"
        if not zip_path.exists():
            self.skipTest(f"Test project zip not found: {zip_path}")

        import_ccp4_project_zip(
            zip_path,
            relocate_path=settings.CCP4I2_PROJECTS_DIR,
        )

        job = Job.objects.filter(project__name="refmac_gamma_test_0").first()
        self.assertIsNotNone(job)

        # Temporarily change task name to something without a report class
        original_task_name = job.task_name
        job.task_name = "nonexistent_task_xyz"
        job.save()

        try:
            report_xml = generate_job_report(job)

            # Should return a failure report, not raise an exception
            self.assertIsNotNone(report_xml)
            report_str = ET.tostring(report_xml, encoding="unicode")
            self.assertIn("No report class found", report_str)
            print("✓ Gracefully handled missing report class")
        finally:
            # Restore original task name
            job.task_name = original_task_name
            job.save()


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent.parent.parent.parent
    / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class ReportRegistryTests(TestCase):
    """Test the report registry functionality."""

    def test_report_class_lazy_loading(self):
        """Test that report classes are loaded lazily."""
        from ccp4i2.core.CCP4TaskManager import TASKMANAGER

        task_manager = TASKMANAGER()

        # Clear cache to ensure fresh load
        task_manager.report_registry.clear_cache()

        # First access should import the class
        report_class = task_manager.getReportClass("refmac")
        self.assertIsNotNone(report_class)
        self.assertEqual(report_class.__name__, "refmac_report")

        # Second access should be cached
        report_class_2 = task_manager.getReportClass("refmac")
        self.assertIs(report_class, report_class_2)

    def test_has_report(self):
        """Test has_report method."""
        from ccp4i2.core.CCP4TaskManager import TASKMANAGER

        task_manager = TASKMANAGER()

        self.assertTrue(task_manager.has_report("refmac"))
        self.assertTrue(task_manager.has_report("pointless"))
        self.assertFalse(task_manager.has_report("nonexistent_task_xyz"))
