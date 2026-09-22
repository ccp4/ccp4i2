"""
Tests for the modern AsyncDatabaseHandler.

These tests demonstrate the clean async API for database-backed job tracking.
"""

import asyncio
import pytest
import uuid
from pathlib import Path
from unittest.mock import Mock, AsyncMock, patch, MagicMock

from ccp4i2.core.CCP4PluginScript import CPluginScript

# NB: this module used to stub Django out at import time --
#     sys.modules['django'] = MagicMock() and friends -- so that the handler
#     could be imported before Django was configured. pytest runs one process,
#     so that replaced Django for the *whole session*: pytest-django's
#     collection hook then did `from django.test import TestCase`, got a
#     MagicMock, and raised "'django' is not a package" as an INTERNALERROR
#     that aborted every run including this directory. `pytest ccp4i2/tests/`
#     collected 3795 tests and then died.
#
#     It is also unnecessary: pytest.ini pins --ds=ccp4i2.config.test_settings,
#     so Django is configured before collection starts. Do not reintroduce it.
#     If a future test needs a module stubbed, scope it with
#     monkeypatch.setitem(sys.modules, ...) so it unwinds after the test.


class TestAsyncDatabaseHandlerDesign:
    """
    Test the design and API of the AsyncDatabaseHandler.

    Note: These are design tests that mock the Django layer.
    Full integration tests require a Django test environment.
    """

    def test_handler_api_surface(self):
        """Test that the handler has the expected modern API"""
        # Import after mocking
        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler

        project_uuid = uuid.uuid4()
        handler = AsyncDatabaseHandler(project_uuid)

        # Check key attributes
        assert handler.project_uuid == project_uuid
        assert handler._project is None  # Lazy loading

        # Check key methods exist with proper signatures
        assert hasattr(handler, 'get_project')
        assert hasattr(handler, 'create_job')
        assert hasattr(handler, 'update_job_status')
        assert hasattr(handler, 'register_output_file')
        assert hasattr(handler, 'register_input_file')
        assert hasattr(handler, 'glean_job_files')
        assert hasattr(handler, 'track_job')

    def test_handler_uses_type_hints(self):
        """Verify that the handler has proper type hints"""
        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler
        import inspect

        # Check constructor
        sig = inspect.signature(AsyncDatabaseHandler.__init__)
        assert sig.parameters['project_uuid'].annotation == uuid.UUID

        # Check create_job
        sig = inspect.signature(AsyncDatabaseHandler.create_job)
        assert sig.parameters['task_name'].annotation == str
        # Optional parameters should have Optional[...] annotation
        assert 'Optional' in str(sig.parameters['title'].annotation)

    def test_plugin_status_conversion(self):
        """Test conversion from CPluginScript status to Job.Status"""
        from ccp4i2.db.async_db_handler import plugin_status_to_job_status
        from ccp4i2.db import models

        # Mock the models.Job.Status enum
        models.Job.Status = Mock()
        models.Job.Status.FINISHED = 6
        models.Job.Status.FAILED = 5
        models.Job.Status.INTERRUPTED = 4
        models.Job.Status.UNSATISFACTORY = 10
        models.Job.Status.TO_DELETE = 9

        # Test conversions
        assert plugin_status_to_job_status(CPluginScript.SUCCEEDED) == 6
        assert plugin_status_to_job_status(CPluginScript.FAILED) == 5
        assert plugin_status_to_job_status(CPluginScript.INTERRUPTED) == 4
        assert plugin_status_to_job_status(CPluginScript.UNSATISFACTORY) == 10
        assert plugin_status_to_job_status(CPluginScript.MARK_TO_DELETE) == 9

        # Unknown status defaults to FAILED
        assert plugin_status_to_job_status(999) == 5

    @pytest.mark.asyncio
    async def test_context_manager_pattern(self):
        """Test that track_job() can be used as an async context manager"""
        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler

        handler = AsyncDatabaseHandler(uuid.uuid4())
        # Not spec=CPluginScript: track_job reads attributes a plugin carries
        # per-instance (name, workDirectory) rather than on the class, which a
        # spec'd Mock refuses. This mirrors the current API -- the old version
        # of this test drove getTaskName() and a statusChanged Qt signal, both
        # of which the Qt-free architecture removed.
        plugin = Mock()
        plugin.get_db_job_id.return_value = None  # so track_job creates the job
        plugin.TASKNAME = "test_plugin"
        plugin.name = "test_job"
        plugin.parent.return_value = None
        plugin.get_status.return_value = CPluginScript.SUCCEEDED
        plugin.container.outputData = None

        # Mock the async methods
        with patch.object(handler, 'create_job', new_callable=AsyncMock) as mock_create:
            with patch.object(handler, 'update_job_status', new_callable=AsyncMock) as mock_update:
                mock_job = Mock()
                mock_job.uuid = uuid.uuid4()
                mock_job.number = "1"
                # track_job does Path(job.directory) to set the plugin's work
                # directory, and Path() will not take a Mock.
                mock_job.directory = "/tmp/track_job_test"
                mock_create.return_value = mock_job

                # Use context manager
                async with handler.track_job(plugin):
                    # Inside context, job should be created and status updated to RUNNING
                    pass

                # Verify create_job was called
                mock_create.assert_called_once()

                # Verify update_job_status was called for RUNNING
                assert mock_update.call_count >= 1

    def test_async_handler_api(self):
        """Verify AsyncDatabaseHandler exposes both modern and legacy APIs."""
        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler
        handler = AsyncDatabaseHandler(uuid.uuid4())

        # Modern snake_case API
        assert hasattr(handler, 'create_job')
        assert hasattr(handler, 'update_job_status')
        assert hasattr(handler, 'track_job')
        assert hasattr(handler, 'register_output_file')
        assert hasattr(handler, 'register_input_file')
        assert hasattr(handler, 'glean_job_files')

        # Legacy camelCase wrappers (used by CPluginScript)
        assert hasattr(handler, 'createSubJob')
        assert hasattr(handler, 'updateJobStatus')
        assert hasattr(handler, 'getProjectDirectory')
        assert hasattr(handler, 'get_file_path_sync')

    def test_job_number_generation_logic(self):
        """
        Test the documented job numbering logic.

        Job numbers follow the pattern:
        - Top-level jobs: "1", "2", "3", ...
        - Child jobs: "1.1", "1.2", "1.3", ...
        - Nested children: "1.1.1", "1.1.2", ...
        """
        # This is a design test - implementation would handle:
        # 1. Parent job with number "1"
        #    - First child: "1.1"
        #    - Second child: "1.2"
        # 2. Child job "1.1"
        #    - First grandchild: "1.1.1"
        # 3. Top-level job increments project.last_job_number

        # The logic is in AsyncDatabaseHandler.create_job()
        # It parses parent.number and appends a new child number
        pass

    def test_status_is_read_from_the_plugin_not_a_signal(self):
        """The Qt-free architecture polls the plugin, it does not connect.

        This test used to assert ``hasattr(CPluginScript, 'statusChanged')``
        and describe track_job connecting to, and disconnecting from, a Qt
        signal. That signal went with PySide2: track_job now reads
        ``plugin.get_status()`` after the ``yield`` and converts it with
        ``plugin_status_to_job_status``. The assertion had been false since
        the migration, but the module could not be collected, so nobody saw
        it fail.
        """
        assert not hasattr(CPluginScript, 'statusChanged'), (
            "a Qt signal has reappeared on CPluginScript; track_job reads "
            "get_status() and does not connect to signals"
        )
        assert hasattr(CPluginScript, 'get_status')


@pytest.mark.integration
class TestAsyncDatabaseHandlerIntegration:
    """
    Integration tests that require a real Django database.

    These tests are marked with @pytest.mark.integration and are skipped
    unless running in a Django test environment.
    """

    @pytest.mark.skip(reason="Requires Django test database setup")
    @pytest.mark.asyncio
    async def test_full_job_lifecycle(self):
        """
        Test complete job lifecycle with real database.

        This test would:
        1. Create a project
        2. Create a job
        3. Track status changes
        4. Register files
        5. Verify database state
        """
        # This requires Django test setup
        pass

    @pytest.mark.skip(reason="Requires Django test database setup")
    @pytest.mark.asyncio
    async def test_nested_job_creation(self):
        """
        Test creating nested jobs with automatic job numbering.
        """
        # This requires Django test setup
        pass


class TestUsagePatterns:
    """
    Tests that demonstrate recommended usage patterns.

    These are example-based tests that show how to use the handler correctly.
    """

    def test_basic_usage_pattern(self):
        """
        Document the basic usage pattern for the new handler.
        """
        example_code = """
        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler

        # Create handler for a project
        handler = AsyncDatabaseHandler(project_uuid=my_project_uuid)

        # Create a plugin
        plugin = CPluginScript(taskName="ctruncate", name="my_job")

        # Track job automatically
        async with handler.track_job(plugin):
            # Job created in database
            # Status updated automatically
            await plugin.execute()
            # Files registered automatically
        """

        # This is documentation - the pattern is clear and simple
        assert "track_job" in example_code
        assert "async with" in example_code

    def test_manual_tracking_pattern(self):
        """
        Document the manual tracking pattern for more control.
        """
        example_code = """
        # For more control, use manual methods
        handler = AsyncDatabaseHandler(project_uuid=my_project_uuid)

        # Create job explicitly
        job = await handler.create_job(
            task_name="ctruncate",
            title="My custom job",
        )

        # Update status manually
        await handler.update_job_status(job.uuid, models.Job.Status.RUNNING)

        # Execute plugin
        result = await plugin.execute()

        # Register files manually
        for output_file in plugin.outputData.files:
            await handler.register_output_file(
                job_uuid=job.uuid,
                file_path=output_file.path,
                file_type=output_file.type,
                param_name=output_file.param_name,
            )

        # Update final status
        if result == CPluginScript.SUCCEEDED:
            await handler.update_job_status(job.uuid, models.Job.Status.FINISHED)
        """

        # Manual pattern gives more control but requires more code
        assert "create_job" in example_code
        assert "update_job_status" in example_code
        assert "register_output_file" in example_code


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
