"""
Tests for CPluginScript base class.
"""

import pytest
import os
from pathlib import Path
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4TaskManager import TASKMANAGER


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
    from ccp4i2.core.base_object.base_classes import CData
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
    from ccp4i2.core.base_object.base_classes import CContainer
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


# ============================================================================
# Tests for dbHandler propagation
# ============================================================================

class TestDbHandlerPropagation:
    """Tests for database handler propagation through plugin hierarchy."""

    def test_dbhandler_captured_from_kwargs(self):
        """Test that dbHandler passed to __init__ is properly captured.

        This is a regression test for a bug where dbHandler was passed via
        **kwargs but never extracted, leaving _dbHandler as None.
        """
        # Create a mock database handler
        class MockDbHandler:
            def __init__(self):
                self.project_uuid = "test-project-uuid"

        mock_handler = MockDbHandler()

        # Create plugin with dbHandler in kwargs
        script = CPluginScript(name="test_script", dbHandler=mock_handler)

        # The handler should be captured
        assert script._dbHandler is not None, \
            "dbHandler passed to __init__ should be captured in _dbHandler"
        assert script._dbHandler is mock_handler, \
            "Captured dbHandler should be the same object that was passed"

    def test_dbhandler_none_when_not_provided(self):
        """Test that _dbHandler is None when not provided."""
        script = CPluginScript(name="test_script")

        assert script._dbHandler is None, \
            "_dbHandler should be None when not provided"

    def test_dbhandler_accessible_from_child_cdatafile(self):
        """Test that CDataFile can access dbHandler via parent hierarchy.

        This tests the full chain: CDataFile -> container -> plugin -> _dbHandler
        """
        from ccp4i2.core.base_object.cdata_file import CDataFile

        class MockDbHandler:
            def __init__(self):
                self.project_uuid = "test-project-uuid"

            def get_file_path_sync(self, file_uuid):
                return f"/mock/path/{file_uuid}"

        mock_handler = MockDbHandler()

        # Create plugin with dbHandler
        script = CPluginScript(name="test_script", dbHandler=mock_handler)

        # Create a CDataFile as a child of inputData
        test_file = CDataFile(parent=script.container.inputData, name="TEST_FILE")

        # The file should be able to find the plugin via _find_plugin_parent()
        found_plugin = test_file._find_plugin_parent()
        assert found_plugin is script, \
            "CDataFile should find the plugin in its parent hierarchy"

        # The file should be able to get the dbHandler via _get_db_handler()
        found_handler = test_file._get_db_handler()
        assert found_handler is mock_handler, \
            "CDataFile._get_db_handler() should return the plugin's _dbHandler"

    def test_setdbdata_also_sets_handler(self):
        """Test that setDbData() can also set the dbHandler."""
        class MockDbHandler:
            pass

        mock_handler = MockDbHandler()

        script = CPluginScript(name="test_script")
        assert script._dbHandler is None

        # setDbData should set the handler
        script.setDbData(handler=mock_handler)

        assert script._dbHandler is mock_handler, \
            "setDbData() should set _dbHandler"

    def test_dbhandler_from_kwargs_not_passed_to_parent(self):
        """Test that dbHandler is popped from kwargs and not passed to parent class.

        If dbHandler is not popped, it would be passed to CData.__init__() via
        super().__init__(**kwargs), which might cause unexpected behavior.
        """
        class MockDbHandler:
            pass

        mock_handler = MockDbHandler()

        # This should not raise any errors about unexpected kwargs
        script = CPluginScript(name="test_script", dbHandler=mock_handler)

        # Verify it was captured
        assert script._dbHandler is mock_handler

    def test_getfullpath_uses_dbhandler_for_dbfileid_resolution(self):
        """Test that CDataFile.getFullPath() uses dbHandler when dbFileId is set.

        This is the key integration test: when a file has dbFileId set (from
        autopopulation), getFullPath() should use the dbHandler to resolve the
        path. Without the dbHandler fix, this would fail and return just the
        baseName instead of the full path.
        """
        from ccp4i2.core.base_object.cdata_file import CDataFile
        import uuid

        test_file_uuid = str(uuid.uuid4())
        expected_path = "/projects/test/CCP4_JOBS/job_1/output.mtz"

        class MockDbHandler:
            def __init__(self):
                self.project_uuid = "test-project-uuid"

            def get_file_path_sync(self, file_uuid):
                """Return the full path for a given file UUID."""
                if str(file_uuid) == test_file_uuid:
                    return expected_path
                return None

        mock_handler = MockDbHandler()

        # Create plugin with dbHandler
        script = CPluginScript(name="test_script", dbHandler=mock_handler)

        # Create a CDataFile as a child of inputData
        test_file = CDataFile(parent=script.container.inputData, name="TEST_FILE")

        # Set file attributes as would happen during autopopulation
        test_file.set({
            "baseName": "output.mtz",
            "relPath": "CCP4_JOBS/job_1",
            "dbFileId": test_file_uuid,
        })

        # getFullPath() should use dbHandler.get_file_path_sync() to resolve
        full_path = test_file.getFullPath()

        assert full_path == expected_path, \
            f"getFullPath() should return '{expected_path}' via dbHandler, got '{full_path}'"

    def test_getfullpath_fails_without_dbhandler(self):
        """Test that path resolution fails gracefully without dbHandler.

        When dbHandler is not available, getFullPath() falls back to baseName,
        which is incorrect for autopopulated files. This test documents the
        expected fallback behavior.
        """
        from ccp4i2.core.base_object.cdata_file import CDataFile
        import uuid

        test_file_uuid = str(uuid.uuid4())

        # Create plugin WITHOUT dbHandler (simulates the bug)
        script = CPluginScript(name="test_script")  # No dbHandler!

        # Verify dbHandler is None
        assert script._dbHandler is None

        # Create a CDataFile as a child of inputData
        test_file = CDataFile(parent=script.container.inputData, name="TEST_FILE")

        # Set file attributes as would happen during autopopulation
        test_file.set({
            "baseName": "output.mtz",
            "relPath": "CCP4_JOBS/job_1",
            "dbFileId": test_file_uuid,
        })

        # Without dbHandler, getFullPath() cannot resolve via dbFileId
        # It will fall back to baseName (which is wrong for autopopulated files)
        full_path = test_file.getFullPath()

        # The path should just be the baseName since dbHandler is not available
        # and the file doesn't actually exist on disk
        assert full_path == "output.mtz", \
            f"Without dbHandler, should fall back to baseName, got '{full_path}'"
