"""
Tests for context-aware parameter setting (CContainer.set_parameter()).

Validates that CContainer.set_parameter() automatically detects CPluginScript
parent context and enables database synchronization when appropriate.
"""

import pytest
from core.base_object.ccontainer import CContainer
from core.base_object.fundamental_types import CInt, CString
from core.base_object.cdata_file import CDataFile


class TestContainerSetParameter:
    """Test CContainer.set_parameter() method."""

    def test_set_parameter_simple(self):
        """Test basic parameter setting without database context."""
        container = CContainer(name="root")
        container.inputData = CContainer(name="inputData")
        container.inputData.NCYCLES = CInt(5, name="NCYCLES")

        # Set parameter using new method
        result = container.set_parameter("inputData.NCYCLES", 10, skip_first=False)

        assert result.value == 10
        assert container.inputData.NCYCLES.value == 10

    def test_set_parameter_with_skip_first(self):
        """Test parameter setting with skip_first for legacy paths."""
        container = CContainer(name="root")
        container.controlParameters = CContainer(name="controlParameters")
        container.controlParameters.NCYCLES = CInt(5, name="NCYCLES")

        # Legacy path includes task name as first element
        result = container.set_parameter("task_name.controlParameters.NCYCLES", 15, skip_first=True)

        assert result.value == 15
        assert container.controlParameters.NCYCLES.value == 15

    def test_set_parameter_file_object(self):
        """Test setting file parameter without database context."""
        container = CContainer(name="root")
        container.inputData = CContainer(name="inputData")
        container.inputData.XYZIN = CDataFile(name="XYZIN")

        # Set file path
        result = container.set_parameter("inputData.XYZIN", "/path/to/file.pdb", skip_first=False)

        # Should set the file path
        assert container.inputData.XYZIN.getFullPath() == "/path/to/file.pdb"

    def test_set_parameter_nested_path(self):
        """Test parameter setting with deeply nested path."""
        container = CContainer(name="root")
        container.level1 = CContainer(name="level1")
        container.level1.level2 = CContainer(name="level2")
        container.level1.level2.value = CString("initial", name="value")

        # Set deeply nested value
        result = container.set_parameter("level1.level2.value", "updated", skip_first=False)

        assert result.value == "updated"
        assert container.level1.level2.value.value == "updated"

    def test_set_parameter_invalid_path(self):
        """Test error handling for invalid path."""
        container = CContainer(name="root")
        container.inputData = CContainer(name="inputData")

        # Try to set non-existent path
        with pytest.raises(AttributeError):
            container.set_parameter("inputData.NONEXISTENT", "value", skip_first=False)


class TestPluginParentDetection:
    """Test _find_plugin_parent() helper method."""

    def test_find_plugin_parent_none(self):
        """Test that standalone container has no plugin parent."""
        container = CContainer(name="standalone")
        plugin_parent = container._find_plugin_parent()

        assert plugin_parent is None

    def test_find_plugin_parent_with_plugin(self):
        """Test finding CPluginScript parent when it exists."""
        # Import CPluginScript
        try:
            from core.CCP4PluginScript import CPluginScript
        except ImportError:
            pytest.skip("CPluginScript not available")

        # Create a plugin with container
        plugin = CPluginScript(name="test_plugin", dummy=True)

        # The container should find its plugin parent
        plugin_parent = plugin.container._find_plugin_parent()

        assert plugin_parent is plugin

    def test_find_plugin_parent_nested_containers(self):
        """Test finding plugin parent through nested container hierarchy."""
        try:
            from core.CCP4PluginScript import CPluginScript
        except ImportError:
            pytest.skip("CPluginScript not available")

        # Create plugin with nested containers
        plugin = CPluginScript(name="test_plugin", dummy=True)

        # Create nested containers with proper parent registration
        input_data = CContainer(name="inputData")
        input_data.set_parent(plugin.container)
        plugin.container.inputData = input_data

        sub_data = CContainer(name="subData")
        sub_data.set_parent(input_data)
        input_data.subData = sub_data

        # Even deeply nested containers should find the plugin
        plugin_parent = sub_data._find_plugin_parent()

        assert plugin_parent is plugin


class TestDatabaseAwareParameterSetting:
    """Test automatic database synchronization when in CPluginScript context."""

    def test_set_parameter_without_dbhandler(self):
        """Test that setting parameter works without dbHandler (no exception)."""
        try:
            from core.CCP4PluginScript import CPluginScript
        except ImportError:
            pytest.skip("CPluginScript not available")

        # Create plugin without dbHandler
        plugin = CPluginScript(name="test_plugin", dummy=True)
        plugin.container.NCYCLES = CInt(5, name="NCYCLES")

        # Should work fine, just won't trigger database update
        result = plugin.container.set_parameter("NCYCLES", 10, skip_first=False)

        assert result.value == 10
        assert plugin.container.NCYCLES.value == 10

    def test_set_parameter_with_mock_dbhandler(self):
        """Test that dbHandler.updateJobStatus() is called when available."""
        try:
            from core.CCP4PluginScript import CPluginScript
        except ImportError:
            pytest.skip("CPluginScript not available")

        # Create plugin with mock dbHandler
        plugin = CPluginScript(name="test_plugin", dummy=True)
        plugin.container.NCYCLES = CInt(5, name="NCYCLES")

        # Mock dbHandler
        class MockDbHandler:
            def __init__(self):
                self.update_called = False
                self.update_job_id = None

            def updateJobStatus(self, jobId, container):
                self.update_called = True
                self.update_job_id = jobId

        mock_handler = MockDbHandler()
        plugin._dbHandler = mock_handler
        plugin._dbJobId = "test-job-123"

        # Set parameter - should trigger database update
        result = plugin.container.set_parameter("NCYCLES", 10, skip_first=False)

        assert result.value == 10
        assert plugin.container.NCYCLES.value == 10
        # Verify database update was triggered
        assert mock_handler.update_called
        assert mock_handler.update_job_id == "test-job-123"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
