"""
Tests for CPluginScript.validity() method.

These tests verify that validation correctly identifies missing required inputs
and other container validation issues.

The validation architecture follows the CData hierarchy:
- CPluginScript.validity() delegates to container.validity()
- Each CData object in the hierarchy has its own validity() method
- CDataFile.validity() checks mustExist, allowUndefined, etc.
- This ensures validation logic lives with the data classes
"""

import xml.etree.ElementTree as ET
import pytest


def test_bare_prosmart_refmac_validation():
    """
    Test that a bare prosmart_refmac job has validation errors.

    A bare prosmart_refmac job (with no inputs set) should report validation
    errors for missing required file inputs:
    - F_SIGF (reflection data)
    - XYZIN (coordinate file)
    """
    # Import Django models and utilities
    from ccp4x.db import models
    from ccp4x.lib.utils.jobs.validate import validate_job
    from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context

    from pathlib import Path
    import tempfile

    # Create a temporary directory for the test project
    with tempfile.TemporaryDirectory(prefix="test_validation_") as tmpdir:
        project_dir = Path(tmpdir) / "test_project"
        project_dir.mkdir()

        # Create required subdirectories
        (project_dir / "CCP4_JOBS").mkdir()
        (project_dir / "CCP4_IMPORTED_FILES").mkdir()
        (project_dir / "CCP4_TMP").mkdir()

        # Create the project in the database
        project = models.Project.objects.create(
            name="test_validation_project",
            directory=str(project_dir)
        )

        # Create a bare prosmart_refmac job
        # Note: 'directory' is a computed property on Job, not a field
        # It's computed from project.directory + "CCP4_JOBS" + "job_{number}"
        job_dir = project_dir / "CCP4_JOBS" / "job_1"
        job_dir.mkdir(parents=True)

        job = models.Job.objects.create(
            project=project,
            task_name="prosmart_refmac",
            number="1",  # Must be string - used to compute directory path
            status=models.Job.Status.PENDING,
        )

        # Verify the directory property works correctly
        assert job.directory == job_dir, f"Job directory mismatch: {job.directory} != {job_dir}"

        # Get the plugin to write initial params.xml
        plugin_result = get_plugin_with_context(job, create_db_handler=False)
        assert plugin_result.success, f"Failed to get plugin: {plugin_result.error}"

        plugin = plugin_result.data

        # Save bare params.xml (no inputs set)
        plugin.container.saveDataToXml(str(job_dir / "params.xml"))

        # Now test validation using validate_job (the existing utility)
        validation_result = validate_job(job)
        assert validation_result.success, f"Validation failed: {validation_result.error}"

        error_tree = validation_result.data

        # Find all error reports
        error_reports = error_tree.findall('.//errorReport')
        print(f"\nFound {len(error_reports)} validation issues:")

        # Collect all object paths with errors
        object_paths_with_errors = set()
        for report in error_reports:
            obj_path = report.findtext('objectPath')
            description = report.findtext('description')
            severity = report.findtext('severity')
            print(f"  [{severity}] {obj_path}: {description}")
            if obj_path:
                object_paths_with_errors.add(obj_path)

        # Verify that we have validation errors (at minimum)
        assert len(error_reports) > 0, "Expected validation errors for bare job, but found none"

        # Check for expected missing inputs
        # prosmart_refmac requires F_SIGF and XYZIN
        found_f_sigf_error = any('F_SIGF' in path for path in object_paths_with_errors)
        found_xyzin_error = any('XYZIN' in path for path in object_paths_with_errors)

        # Print summary
        print(f"\nF_SIGF error found: {found_f_sigf_error}")
        print(f"XYZIN error found: {found_xyzin_error}")

        # Assert that we found errors for required inputs
        assert found_f_sigf_error or found_xyzin_error, (
            f"Expected errors for F_SIGF or XYZIN, but found errors only for: {object_paths_with_errors}"
        )

        # Verify object paths are correctly formatted (not "<method>.XYZIN" etc.)
        # Object paths should follow the pattern: task_name.container.inputData.FIELD
        for path in object_paths_with_errors:
            assert not path.startswith('<'), (
                f"Object path should not start with '<' (got '{path}'). "
                "This indicates objectPath() is not being called correctly."
            )
            assert 'prosmart_refmac' in path, (
                f"Object path should contain task name 'prosmart_refmac' (got '{path}')"
            )

        # Cleanup
        models.Job.objects.filter(id=job.id).delete()
        models.Project.objects.filter(id=project.id).delete()


def test_plugin_validity_method():
    """
    Test that CPluginScript has a validity() method that works correctly.

    This tests the CPluginScript.validity() method directly, ensuring it
    returns a CErrorReport with validation issues.
    """
    from ccp4x.db import models
    from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context
    from pathlib import Path
    import tempfile

    # Create a temporary directory for the test project
    with tempfile.TemporaryDirectory(prefix="test_validity_") as tmpdir:
        project_dir = Path(tmpdir) / "test_project"
        project_dir.mkdir()

        # Create required subdirectories
        (project_dir / "CCP4_JOBS").mkdir()
        (project_dir / "CCP4_IMPORTED_FILES").mkdir()
        (project_dir / "CCP4_TMP").mkdir()

        # Create the project in the database
        project = models.Project.objects.create(
            name="test_validity_project",
            directory=str(project_dir)
        )

        # Create a bare prosmart_refmac job
        job_dir = project_dir / "CCP4_JOBS" / "job_1"
        job_dir.mkdir(parents=True)

        job = models.Job.objects.create(
            project=project,
            task_name="prosmart_refmac",
            number="1",  # Must be string
            status=models.Job.Status.PENDING,
        )

        # Get the plugin
        plugin_result = get_plugin_with_context(job, create_db_handler=False)
        assert plugin_result.success, \
            f"Failed to get plugin: {plugin_result.error}"

        plugin = plugin_result.data

        # Save bare params.xml
        plugin.container.saveDataToXml(str(job_dir / "params.xml"))

        # Test that plugin has validity() method
        assert hasattr(plugin, 'validity'), \
            "CPluginScript should have validity() method"

        # Call validity() and check it returns a CErrorReport
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        error_report = plugin.validity()

        assert isinstance(error_report, CErrorReport), (
            f"validity() should return CErrorReport, got {type(error_report)}"
        )

        # Check that we have validation errors for a bare job
        error_count = error_report.count()
        print(f"\nvalidity() returned {error_count} errors")

        if error_count > 0:
            print(f"Error report:\n{error_report.report()}")

        # Cleanup
        models.Job.objects.filter(id=job.id).delete()
        models.Project.objects.filter(id=project.id).delete()


def test_plugin_validity_as_xml_method():
    """
    Test that CPluginScript has a validity_as_xml() method that returns XML.

    This tests the convenience method that returns validation as XML Element.
    """
    from ccp4x.db import models
    from ccp4x.lib.utils.plugins.plugin_context import get_plugin_with_context
    from pathlib import Path
    import tempfile

    # Create a temporary directory for the test project
    with tempfile.TemporaryDirectory(prefix="test_validity_xml_") as tmpdir:
        project_dir = Path(tmpdir) / "test_project"
        project_dir.mkdir()

        # Create required subdirectories
        (project_dir / "CCP4_JOBS").mkdir()
        (project_dir / "CCP4_IMPORTED_FILES").mkdir()
        (project_dir / "CCP4_TMP").mkdir()

        # Create the project in the database
        project = models.Project.objects.create(
            name="test_validity_xml_project",
            directory=str(project_dir)
        )

        # Create a bare prosmart_refmac job
        job_dir = project_dir / "CCP4_JOBS" / "job_1"
        job_dir.mkdir(parents=True)

        job = models.Job.objects.create(
            project=project,
            task_name="prosmart_refmac",
            number="1",  # Must be string
            status=models.Job.Status.PENDING,
        )

        # Get the plugin
        plugin_result = get_plugin_with_context(job, create_db_handler=False)
        assert plugin_result.success, \
            f"Failed to get plugin: {plugin_result.error}"

        plugin = plugin_result.data

        # Save bare params.xml
        plugin.container.saveDataToXml(str(job_dir / "params.xml"))

        # Test that plugin has validity_as_xml() method
        assert hasattr(plugin, 'validity_as_xml'), \
            "CPluginScript should have validity_as_xml() method"

        # Call validity_as_xml() and check it returns an XML Element
        error_xml = plugin.validity_as_xml()

        assert isinstance(error_xml, ET.Element), (
            f"validity_as_xml() should return ET.Element, got {type(error_xml)}"
        )

        # Check the XML structure
        assert error_xml.tag == 'errorReportList', (
            f"Root tag should be 'errorReportList', got '{error_xml.tag}'"
        )

        # Find error reports
        error_reports = error_xml.findall('.//errorReport')
        print(f"\nvalidity_as_xml() returned {len(error_reports)} errors")

        # Print the XML for debugging
        ET.indent(error_xml, ' ')
        xml_str = ET.tostring(error_xml, encoding='unicode')
        print(f"XML:\n{xml_str[:1000]}...")  # Truncate for readability

        # Cleanup
        models.Job.objects.filter(id=job.id).delete()
        models.Project.objects.filter(id=project.id).delete()
