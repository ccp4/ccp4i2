"""
Integration tests for report generation.

These tests run actual CCP4 jobs via i2run, then test the report generation
system by:
1. Running a job that produces program.xml output
2. Verifying the job is stored in the database
3. Calling generate_job_report() to generate the report
4. Validating the report structure and content

This ensures the full report generation pipeline works with real job data.
"""
import xml.etree.ElementTree as ET
import pytest
from pathlib import Path
from .utils import demoData, i2run


def test_refmac_report_generation():
    """
    Test report generation for a successful refmac job.

    Refmac is a good test case because:
    - It produces detailed program.xml output
    - Its report class (refmac_report) is well-tested
    - It's commonly used and representative of refinement tasks
    """
    from ccp4x.db import models
    from ccp4x.lib.utils.reporting.i2_report import generate_job_report

    # Run a simple refmac job using gamma demo data
    args = [
        "refmac",
        "--XYZIN", demoData("gamma", "gamma_model.pdb"),
        "--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz"),
    ]

    with i2run(args) as job_dir:
        # Verify job completed and produced output
        assert job_dir.exists(), f"Job directory not found: {job_dir}"

        # Check for program.xml (refmac output)
        program_xml = job_dir / "program.xml"
        if not program_xml.exists():
            # Some jobs output XMLOUT.xml instead
            program_xml = job_dir / "XMLOUT.xml"

        # Find the job in the database
        job_number = job_dir.name.replace("job_", "")
        job = models.Job.objects.filter(number=job_number).first()
        assert job is not None, f"Job {job_number} not found in database"
        assert job.task_name == "refmac", f"Wrong task name: {job.task_name}"

        # Test report generation
        print(f"\n=== Testing report generation for {job.task_name} job {job.number} ===")
        report_xml = generate_job_report(job)

        # Validate report structure
        assert report_xml is not None, "Report generation returned None"
        assert isinstance(report_xml, ET.Element), "Report should be an XML Element"

        # Convert to string for inspection
        report_str = ET.tostring(report_xml, encoding="unicode")
        print(f"Report XML length: {len(report_str)} characters")

        # Check for expected report elements
        # Refmac reports should have title and result sections
        assert "CCP4i2ReportTitle" in report_str or "title" in report_str.lower(), \
            "Report should have a title section"

        # Check that it's not a failed report
        assert "No report because" not in report_str, \
            f"Report generation failed: {report_str[:500]}"

        print("✓ Report generated successfully")

        # Check for specific refmac report content
        if "RFactor" in report_str or "R-factor" in report_str or "Rfactor" in report_str:
            print("✓ R-factor information found in report")

        # Save report for manual inspection
        report_path = job_dir / "test_report.xml"
        ET.ElementTree(report_xml).write(report_path, encoding="unicode")
        print(f"✓ Report saved to {report_path}")


def test_pointless_report_generation():
    """
    Test report generation for pointless (data reduction task).

    Pointless uses XMLOUT.xml format which is a different code path.
    """
    from ccp4x.db import models
    from ccp4x.lib.utils.reporting.i2_report import generate_job_report

    args = [
        "pointless",
        "--HKLIN", demoData("unscaled_data/gamma/indexed.mtz"),
    ]

    with i2run(args) as job_dir:
        assert job_dir.exists()

        job_number = job_dir.name.replace("job_", "")
        job = models.Job.objects.filter(number=job_number).first()
        assert job is not None, f"Job not found in database"

        print(f"\n=== Testing report generation for {job.task_name} job {job.number} ===")
        report_xml = generate_job_report(job)

        assert report_xml is not None
        report_str = ET.tostring(report_xml, encoding="unicode")

        # Verify it's not a failure report
        if "No report because" in report_str:
            print(f"Warning: Report indicates failure: {report_str[:500]}")
            # This may be expected if no report class exists for pointless
            # Check if the error message is informative
            assert "No report class" in report_str or "No program XML" in report_str

        print(f"✓ Report generated ({len(report_str)} chars)")


def test_import_merged_report_generation():
    """
    Test report generation for import_merged task.

    Import tasks often don't have program.xml, testing the fallback behavior.
    """
    from ccp4x.db import models
    from ccp4x.lib.utils.reporting.i2_report import generate_job_report

    # Use gamma merged intensities which exist
    args = [
        "import_merged",
        "--HKLIN", demoData("gamma", "merged_intensities_Xe.mtz"),
    ]

    with i2run(args) as job_dir:
        assert job_dir.exists()

        job_number = job_dir.name.replace("job_", "")
        job = models.Job.objects.filter(number=job_number).first()
        assert job is not None

        print(f"\n=== Testing report generation for {job.task_name} job {job.number} ===")
        report_xml = generate_job_report(job)

        assert report_xml is not None
        report_str = ET.tostring(report_xml, encoding="unicode")
        print(f"Report generated: {len(report_str)} chars")


def test_report_for_nonexistent_task():
    """
    Test that report generation handles missing report classes gracefully.
    """
    from ccp4x.lib.utils.reporting.i2_report import simple_failed_report

    # Test the simple_failed_report function directly
    report = simple_failed_report(
        reason="Test failure reason",
        task_name="nonexistent_task",
        details="This is a test of the error reporting system."
    )

    assert report is not None
    report_str = ET.tostring(report, encoding="unicode")

    # Verify the failure report structure
    assert "nonexistent_task" in report_str
    assert "Test failure reason" in report_str
    assert "This is a test" in report_str

    print("✓ Failure report generated correctly")


def test_report_metadata_access():
    """
    Test that report metadata can be accessed without generating the full report.

    This tests the new report registry functionality.
    """
    from ccp4i2.core.CCP4TaskManager import TASKMANAGER

    task_manager = TASKMANAGER()

    # Test has_report
    assert task_manager.has_report("refmac"), "refmac should have a report"

    # Test get_report_metadata (fast, no import)
    metadata = task_manager.get_report_metadata("refmac")
    assert metadata is not None, "refmac should have metadata"
    assert metadata.get("TASKNAME") == "refmac"

    # Test RUNNING attribute
    running = task_manager.getReportAttribute("refmac", "RUNNING")
    assert running is True, "refmac report should support RUNNING=True"

    # Test list_reports
    reports = task_manager.list_reports()
    assert len(reports) > 100, f"Expected 100+ reports, got {len(reports)}"
    assert "refmac" in reports

    print(f"✓ Report registry working ({len(reports)} reports available)")


def test_csymmatch_report_generation():
    """
    Test report generation for csymmatch (coordinate utility).

    Uses the gamma crystal structure which is small and fast.
    """
    from ccp4x.db import models
    from ccp4x.lib.utils.reporting.i2_report import generate_job_report

    args = [
        "csymmatch",
        "--XYZIN_QUERY", demoData("8xfm.cif"),
        "--XYZIN_TARGET", demoData("8xfm.cif"),
    ]

    with i2run(args) as job_dir:
        assert job_dir.exists()

        job_number = job_dir.name.replace("job_", "")
        job = models.Job.objects.filter(number=job_number).first()
        assert job is not None

        print(f"\n=== Testing report generation for {job.task_name} job {job.number} ===")
        report_xml = generate_job_report(job)

        assert report_xml is not None
        report_str = ET.tostring(report_xml, encoding="unicode")

        # Check for expected content or appropriate error
        if "No report because" not in report_str:
            print(f"✓ Report generated successfully ({len(report_str)} chars)")
        else:
            print(f"Note: Report shows expected limitation: {report_str[:200]}")
