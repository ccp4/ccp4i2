"""
Tests for CPluginScript base class.
"""

import pytest
import os
from pathlib import Path
from core.CCP4PluginScript import CPluginScript
from core.CCP4TaskManager import TASKMANAGER


def test_cpluginscript_instantiation():
    """Test basic CPluginScript instantiation."""
    script = CPluginScript(name="test_script")

    assert script.name == "test_script"
    assert script.container is not None
    assert script.container.inputData is not None
    assert script.container.outputData is not None
    assert script.container.controlParameters is not None
    assert script.container.guiAdmin is not None


def test_cpluginscript_subcontainers():
    """Test that sub-containers are properly initialized."""
    script = CPluginScript(name="test_script")

    # Check that sub-containers have correct names (using name attribute)
    assert script.container.inputData.name == "inputData"
    assert script.container.outputData.name == "outputData"
    assert script.container.controlParameters.name == "controlParameters"
    assert script.container.guiAdmin.name == "guiAdmin"

    # Check that they are children of the main container
    assert script.container.inputData.parent == script.container
    assert script.container.outputData.parent == script.container


def test_cpluginscript_error_report():
    """Test error reporting functionality."""
    script = CPluginScript(name="test_script")

    assert script.errorReport is not None
    # Error report is initially empty (no errors)
    assert not script.errorReport  # Empty error report is falsy


def test_cpluginscript_subclass():
    """Test creating a subclass of CPluginScript."""

    class TestWrapper(CPluginScript):
        TASKMODULE = 'utility'
        TASKTITLE = 'Test Wrapper'
        TASKNAME = 'test_wrapper'
        TASKCOMMAND = 'test_command'
        TASKVERSION = 1.0

    wrapper = TestWrapper()

    assert wrapper.TASKMODULE == 'utility'
    assert wrapper.TASKTITLE == 'Test Wrapper'
    assert wrapper.TASKNAME == 'test_wrapper'
    assert wrapper.TASKCOMMAND == 'test_command'
    assert wrapper.TASKVERSION == 1.0
    assert wrapper.name == 'test_wrapper'


def test_cpluginscript_process_workflow():
    """Test the basic process workflow methods exist."""
    script = CPluginScript(name="test_script")

    # Check that all required workflow methods exist
    assert hasattr(script, 'process')
    assert hasattr(script, 'checkInputData')
    assert hasattr(script, 'checkOutputData')
    assert hasattr(script, 'processInputFiles')
    assert hasattr(script, 'makeCommandAndScript')
    assert hasattr(script, 'startProcess')
    assert hasattr(script, 'postProcess')
    assert hasattr(script, 'postProcessCheck')
    assert hasattr(script, 'processOutputFiles')
    assert hasattr(script, 'reportStatus')


def test_cpluginscript_hierarchy():
    """Test that CPluginScript supports the hierarchy system."""
    script = CPluginScript(name="test_script")

    # CPluginScript should inherit from CData
    from core.base_object.base_classes import CData
    assert isinstance(script, CData)

    # Container should be a child of the script
    assert script.container.parent == script

    # Sub-containers should be children of the main container
    assert script.container.inputData.parent == script.container
    assert script.container.outputData.parent == script.container
    assert script.container.controlParameters.parent == script.container
    assert script.container.guiAdmin.parent == script.container

    # Should support the object_path API
    assert script.object_path() == "test_script"
    assert script.container.object_path() == "test_script.container"
    assert script.container.inputData.object_path() == "test_script.container.inputData"


def test_cpluginscript_events():
    """Test that CPluginScript supports the event system."""
    script = CPluginScript(name="test_script")

    # Should have event emitters from HierarchicalObject
    assert hasattr(script, 'child_added')
    assert hasattr(script, 'child_removed')

    # Test event emission
    event_fired = []

    def on_child_added(child):
        event_fired.append(child)

    script.child_added.connect(on_child_added)

    # Adding a new container should fire the event
    from core.base_object.base_classes import CContainer
    new_container = CContainer(parent=script, name="test_container")

    assert len(event_fired) == 1
    assert event_fired[0] == new_container


# ============================================================================
# Tests for .def.xml loading functionality
# ============================================================================

@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
class TestDefXmlLoading:
    """Tests for automatic .def.xml loading via CTaskManager."""

    def test_taskmanager_locate_def_xml(self):
        """Test CTaskManager.locate_def_xml() finds plugin definition files."""
        tm = TASKMANAGER()

        # Test finding a common plugin
        path = tm.locate_def_xml('pointless')

        assert path is not None, "Should find pointless.def.xml"
        assert path.exists(), f"Path should exist: {path}"
        assert path.name == 'pointless.def.xml'
        assert path.is_absolute()

    def test_taskmanager_locate_def_xml_multiple_plugins(self):
        """Test locating multiple different plugins."""
        tm = TASKMANAGER()

        test_plugins = ['refmac', 'pointless', 'aimless', 'coot1']

        for plugin_name in test_plugins:
            path = tm.locate_def_xml(plugin_name)
            assert path is not None, f"Should find {plugin_name}"
            assert path.exists(), f"{plugin_name} path should exist"
            assert plugin_name in path.name

    def test_taskmanager_locate_def_xml_with_version(self):
        """Test locating plugin with specific version."""
        tm = TASKMANAGER()

        # buccaneer_mr has version 0.0 in the lookup
        path = tm.locate_def_xml('buccaneer_mr', version='0.0')

        assert path is not None, "Should find versioned plugin"
        assert path.exists()
        assert path.name == 'buccaneer_mr.def.xml'

    def test_taskmanager_locate_def_xml_not_found(self):
        """Test that non-existent plugins return None."""
        tm = TASKMANAGER()

        path = tm.locate_def_xml('nonexistent_plugin_xyz')

        assert path is None, "Should return None for non-existent plugin"

    def test_pluginscript_loads_own_def_xml(self):
        """Test that CPluginScript subclass automatically loads its .def.xml."""

        class PointlessWrapper(CPluginScript):
            TASKNAME = 'pointless'
            TASKVERSION = None

        plugin = PointlessWrapper()

        # Check that def file was found and loaded
        assert plugin.defFile is not None, "Should have located .def.xml"
        assert Path(plugin.defFile).exists(), "Def file should exist"
        assert 'pointless.def.xml' in plugin.defFile

        # Check that no errors occurred during loading
        assert plugin.errorReport.count() == 0, \
            f"Should load without errors, got: {plugin.errorReport.report()}"

    def test_pluginscript_loads_aimless_def_xml(self):
        """Test loading aimless plugin definition."""

        class AimlessWrapper(CPluginScript):
            TASKNAME = 'aimless'
            TASKVERSION = None

        plugin = AimlessWrapper()

        # Check that def file was loaded
        assert plugin.defFile is not None
        assert 'aimless.def.xml' in plugin.defFile

        # The def file should have populated the containers
        # (Note: inputData.items() may not be available depending on
        # how DefXmlParser populates containers)
        assert plugin.container.inputData is not None
        assert plugin.container.outputData is not None
        assert plugin.container.controlParameters is not None

    def test_pluginscript_loads_refmac_def_xml(self):
        """Test loading refmac plugin definition."""

        class RefmacWrapper(CPluginScript):
            TASKNAME = 'refmac'
            TASKVERSION = None

        plugin = RefmacWrapper()

        # Check that def file was loaded
        assert plugin.defFile is not None
        assert 'refmac.def.xml' in plugin.defFile
        assert Path(plugin.defFile).exists()

        # Should have standard containers
        assert plugin.container is not None
        assert plugin.container.inputData is not None
        assert plugin.container.outputData is not None

    def test_pluginscript_without_taskname_no_def_load(self):
        """Test that plugin without TASKNAME doesn't try to load .def.xml."""

        class GenericWrapper(CPluginScript):
            TASKNAME = None  # No task name

        plugin = GenericWrapper()

        # Should not have loaded a def file
        assert plugin.defFile is None

        # Should still have basic structure
        assert plugin.container is not None
        assert plugin.container.inputData is not None

    def test_pluginscript_versioned_loading(self):
        """Test loading plugin with specific version."""

        class BuccaneerWrapper(CPluginScript):
            TASKNAME = 'buccaneer_mr'
            TASKVERSION = '0.0'

        plugin = BuccaneerWrapper()

        # Should find the versioned def file
        assert plugin.defFile is not None
        assert 'buccaneer_mr.def.xml' in plugin.defFile

    def test_pluginscript_hierarchy_after_def_load(self):
        """Test that hierarchy is preserved after loading .def.xml."""

        class PointlessWrapper(CPluginScript):
            TASKNAME = 'pointless'
            TASKVERSION = None

        plugin = PointlessWrapper()

        # Check hierarchy is intact
        assert plugin.container.parent == plugin
        assert plugin.container.inputData.parent == plugin.container
        assert plugin.container.outputData.parent == plugin.container

        # Check object paths work
        assert plugin.object_path() == 'pointless'
        assert 'pointless' in plugin.container.object_path()
