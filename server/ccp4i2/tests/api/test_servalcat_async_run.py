"""
Test servalcat_pipe execution via async_run_job code path.

This test exercises the HTTP API code path (run_job_async) for servalcat_pipe,
which is different from the i2run command-line path. The key difference is
that run_job_async calls validity() before process(), allowing plugins to
adjust qualifiers for embedded wrappers.

The servalcat_pipe has embedded metalCoordWrapper which requires XYZIN to be
set when run standalone. In the pipeline context, XYZIN is programmatically
filled, so servalcat_pipe.validity() sets allowUndefined=True.
"""

import logging
from pathlib import Path
from shutil import rmtree
import tempfile

import pytest

from ...db import models

logger = logging.getLogger(f"ccp4i2::{__name__}")


# Create a unique test directory
TEST_DIR = Path(tempfile.gettempdir()) / "ccp4i2_test_servalcat_async"


@pytest.mark.django_db(transaction=True)
class TestServalcatAsyncRun:
    """Test servalcat_pipe execution via async_run_job code path."""

    @pytest.fixture(autouse=True)
    def setup_project(self, db):
        """Set up a test project for each test."""
        import uuid

        # Create test directory
        TEST_DIR.mkdir(parents=True, exist_ok=True)

        # Create project in database
        self.project = models.Project.objects.create(
            name=f"servalcat_async_test_{uuid.uuid4().hex[:8]}",
            directory=str(TEST_DIR / "project"),
        )
        Path(self.project.directory).mkdir(parents=True, exist_ok=True)

        yield

        # Cleanup
        try:
            rmtree(TEST_DIR, ignore_errors=True)
        except Exception:
            pass

    @pytest.mark.asyncio
    async def test_metalCoord_allowUndefined_is_set_before_validation(self, setup_project):
        """
        Unit test to verify that servalcat_pipe.validity() sets allowUndefined=True
        on metalCoordWrapper.inputData.XYZIN BEFORE parent validation runs.

        This is a more focused test that just checks the validity() behavior
        without actually running the job.
        """
        from ccp4i2.core.CCP4TaskManager import CTaskManager

        task_manager = CTaskManager()
        plugin_class = task_manager.get_plugin_class("servalcat_pipe")
        plugin = plugin_class(
            workDirectory=str(TEST_DIR / "test_validity_unit"),
            name="test_validity_unit",
        )

        Path(plugin.workDirectory).mkdir(parents=True, exist_ok=True)

        # Check initial state of metalCoordWrapper.inputData.XYZIN
        metal_coord_xyzin = plugin.container.metalCoordWrapper.inputData.XYZIN
        initial_allowUndefined = metal_coord_xyzin.get_qualifier('allowUndefined')
        logger.info(f"Initial allowUndefined for metalCoordWrapper.inputData.XYZIN: {initial_allowUndefined}")

        # From metalCoord.def.xml, this should be False
        assert initial_allowUndefined is False, \
            f"Expected initial allowUndefined to be False, got {initial_allowUndefined}"

        # Call validity() - this should set allowUndefined=True
        error_report = plugin.validity()

        # Check that allowUndefined was set to True
        after_allowUndefined = metal_coord_xyzin.get_qualifier('allowUndefined')
        logger.info(f"After validity() allowUndefined for metalCoordWrapper.inputData.XYZIN: {after_allowUndefined}")

        assert after_allowUndefined is True, \
            f"Expected allowUndefined to be True after validity(), got {after_allowUndefined}"

        # The error report should NOT contain ERROR for metalCoordWrapper.inputData.XYZIN
        # (it may contain WARNING though, which is fine)
        from ccp4i2.core.base_object.error_reporting import SEVERITY_ERROR
        errors_for_xyzin = [
            e for e in error_report.getErrors()
            if 'metalCoordWrapper.inputData.XYZIN' in e.get('name', '')
               and e.get('severity', 0) >= SEVERITY_ERROR
        ]

        assert len(errors_for_xyzin) == 0, \
            f"Should NOT have ERROR for metalCoordWrapper.inputData.XYZIN after validity(), but found: {errors_for_xyzin}"

        logger.info("SUCCESS: validity() correctly sets allowUndefined=True before validation")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
