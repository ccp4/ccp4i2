"""Tests for CPluginScript.makePluginObject method."""

import pytest
import os
from ccp4i2.core.CCP4PluginScript import CPluginScript


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
class TestMakePluginObject:
    """Tests for the makePluginObject method."""

    def test_makepluginobject_creates_instance(self):
        """Test that makePluginObject creates a valid plugin instance."""
        # Create a parent plugin
        parent = CPluginScript(name="parent_task")

        # Create a sub-plugin using makePluginObject
        sub_plugin = parent.makePluginObject("pointless")

        assert sub_plugin is not None, "Should create plugin instance"
        assert isinstance(sub_plugin, CPluginScript), "Should be a CPluginScript instance"
        assert sub_plugin.TASKNAME == "pointless", "Should have correct TASKNAME"

    def test_makepluginobject_with_version(self):
        """Test that makePluginObject respects version parameter."""
        parent = CPluginScript(name="parent_task")

        # Create sub-plugin with specific version
        sub_plugin = parent.makePluginObject("pointless", version=0.0)

        assert sub_plugin is not None, "Should create plugin instance with version"
        assert sub_plugin.TASKNAME == "pointless"

    def test_makepluginobject_nonexistent_plugin(self):
        """Test that makePluginObject handles non-existent plugins."""
        parent = CPluginScript(name="parent_task")

        # Try to create non-existent plugin
        sub_plugin = parent.makePluginObject("nonexistent_task_xyz")

        assert sub_plugin is None, "Should return None for non-existent plugin"
        # Check that error was logged (errorReport should have errors)
        assert bool(parent.errorReport), "Should log error"

    def test_makepluginobject_passes_kwargs(self):
        """Test that makePluginObject passes kwargs to the plugin constructor."""
        parent = CPluginScript(name="parent_task")

        # Create sub-plugin with custom name
        sub_plugin = parent.makePluginObject("pointless", name="custom_pointless")

        assert sub_plugin is not None, "Should create plugin instance"
        assert sub_plugin.name == "custom_pointless", "Should use custom name"

    def test_makepluginobject_multiple_plugins(self):
        """Test creating multiple different plugins."""
        parent = CPluginScript(name="parent_task")

        # Create multiple sub-plugins (using plugins with minimal dependencies)
        plugin1 = parent.makePluginObject("pointless")
        plugin2 = parent.makePluginObject("aimless")
        plugin3 = parent.makePluginObject("unique")

        assert plugin1 is not None and plugin1.TASKNAME == "pointless"
        assert plugin2 is not None and plugin2.TASKNAME == "aimless"
        # unique plugin may not exist - just check it doesn't crash
        # The important thing is that we can call makePluginObject multiple times

    def test_makepluginobject_containers_initialized(self):
        """Test that created plugins have proper container structure."""
        parent = CPluginScript(name="parent_task")

        sub_plugin = parent.makePluginObject("pointless")

        assert sub_plugin is not None
        assert hasattr(sub_plugin, 'container'), "Should have container"
        # Check that container has the expected sub-containers
        assert hasattr(sub_plugin.container, 'inputData'), "Should have container.inputData"
        assert hasattr(sub_plugin.container, 'outputData'), "Should have container.outputData"
        assert hasattr(sub_plugin.container, 'controlParameters'), "Should have container.controlParameters"
