"""
Test CPluginScript async execution with signal system.

This verifies that our CPluginScript now works exactly like the Qt version
for async plugin chaining (the demo_copycell pattern).
"""

import pytest
import time
from pathlib import Path
import tempfile
from core.CCP4PluginScript import CPluginScript
from core.base_object.signal_system import Slot


class SimpleAsyncPlugin(CPluginScript):
    """
    Simple async plugin for testing.

    Mimics mtzdump/pdbset - runs a shell command asynchronously.
    """

    TASKMODULE = 'test'
    TASKTITLE = 'Simple Async Test Plugin'
    TASKNAME = 'simple_async'
    TASKCOMMAND = 'echo'
    TASKVERSION = '1.0'
    ASYNCHRONOUS = True  # Enable async execution

    def __init__(self, message: str = "test output", **kwargs):
        super().__init__(**kwargs)
        self.message = message

    def process(self):
        """Start the async process."""
        # Set command line
        self.commandLine = [self.message]

        # Call parent's process() which will use async execution
        return super().process()


class SimplePipeline(CPluginScript):
    """
    Simple pipeline that chains two async plugins.

    This mimics the demo_copycell pattern:
    1. Run plugin1 (like mtzdump)
    2. When it finishes, run plugin2 (like pdbset)
    3. Report final status
    """

    TASKMODULE = 'test'
    TASKTITLE = 'Simple Pipeline Test'
    TASKNAME = 'simple_pipeline'
    TASKVERSION = '1.0'
    ASYNCHRONOUS = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def process(self):
        """Start the pipeline."""
        # Create first plugin
        self.plugin1 = SimpleAsyncPlugin(
            message="plugin1 output",
            name="plugin1",
            workDirectory=self.workDirectory
        )

        # Connect finished signal to handler
        self.connectSignal(self.plugin1, 'finished', self.on_plugin1_finished)

        # Start plugin1
        self.plugin1.process()

        return self.SUCCEEDED

    @Slot(dict)
    def on_plugin1_finished(self, status_dict):
        """Called when plugin1 finishes."""
        print(f"Pipeline: plugin1 finished with status {status_dict['finishStatus']}")

        status = status_dict['finishStatus']
        if status == self.FAILED:
            self.reportStatus(self.FAILED)
            return

        # Create second plugin
        self.plugin2 = SimpleAsyncPlugin(
            message="plugin2 output",
            name="plugin2",
            workDirectory=self.workDirectory
        )

        # Connect to final handler
        self.connectSignal(self.plugin2, 'finished', self.on_pipeline_finished)

        # Start plugin2
        self.plugin2.process()

    @Slot(dict)
    def on_pipeline_finished(self, status_dict):
        """Called when plugin2 finishes."""
        print(f"Pipeline: plugin2 finished with status {status_dict['finishStatus']}")

        # Report final status
        self.reportStatus(status_dict['finishStatus'])


class TestCPluginScriptAsync:
    """Test CPluginScript async execution."""

    def test_single_async_plugin(self):
        """Test single plugin with async execution."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = {'status': None, 'completed': False}

            @Slot(dict)
            def on_finished(status_dict):
                result['status'] = status_dict['finishStatus']
                result['completed'] = True

            # Create plugin
            plugin = SimpleAsyncPlugin(
                message="async test",
                workDirectory=tmpdir,
                name="test_async"
            )

            # Connect signal
            plugin.finished.connect(on_finished, weak=False)

            # Start plugin
            status = plugin.process()
            assert status == CPluginScript.SUCCEEDED

            # Wait for completion
            for _ in range(100):  # Max 10 seconds
                if result['completed']:
                    break
                time.sleep(0.1)

            # Verify
            assert result['completed'], "Plugin did not complete"
            assert result['status'] == CPluginScript.SUCCEEDED

    def test_pipeline_two_async_plugins(self):
        """Test pipeline with two chained async plugins (demo_copycell pattern)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = {'status': None, 'completed': False}

            @Slot(dict)
            def on_pipeline_finished(status_dict):
                result['status'] = status_dict['finishStatus']
                result['completed'] = True

            # Create pipeline
            pipeline = SimplePipeline(
                workDirectory=tmpdir,
                name="test_pipeline"
            )

            # Connect to final signal
            pipeline.finished.connect(on_pipeline_finished, weak=False)

            # Start pipeline
            status = pipeline.process()
            assert status == CPluginScript.SUCCEEDED

            # Wait for completion (both plugins need to finish)
            for _ in range(100):  # Max 10 seconds
                if result['completed']:
                    break
                time.sleep(0.1)

            # Verify
            assert result['completed'], "Pipeline did not complete"
            assert result['status'] == CPluginScript.SUCCEEDED

            # Verify both plugins were created
            assert hasattr(pipeline, 'plugin1'), "plugin1 was not created"
            assert hasattr(pipeline, 'plugin2'), "plugin2 was not created"

    def test_connect_signal_api(self):
        """Test that connectSignal() API works correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            plugin = SimpleAsyncPlugin(
                message="signal test",
                workDirectory=tmpdir
            )

            # Test connectSignal
            called = {'value': False}

            def handler(status_dict):
                called['value'] = True

            # This should work without error
            plugin.connectSignal(plugin, 'finished', handler)

            # Manually emit signal to test connection
            plugin.finished.emit({'finishStatus': CPluginScript.SUCCEEDED})

            assert called['value'], "Signal handler was not called"

    def test_async_flag_detection(self):
        """Test that ASYNCHRONOUS flag is detected."""
        plugin = SimpleAsyncPlugin()
        assert plugin.ASYNCHRONOUS == True

        # Verify _startProcessAsync would be called
        # (We can't easily test the actual routing without running a process)


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
