"""
Integration tests for async plugin execution with shell commands.

This test bootstraps from simple shell commands to verify the async
execution infrastructure works before testing with real plugins.
"""

import pytest
import asyncio
import time
import tempfile
from pathlib import Path
from ccp4i2.core.base_object.signal_system import Signal, Slot
from ccp4i2.core.async_process_manager import AsyncProcessManager


class MockPlugin:
    """
    Mock plugin that mimics CPluginScript async behavior.

    This is a simplified version that tests the core async pattern:
    - Has a 'finished' signal
    - Calls process() which starts async subprocess
    - Emits 'finished' signal when subprocess completes
    """

    # Status codes (like CPluginScript)
    SUCCEEDED = 0
    FAILED = 1

    def __init__(self, name: str, command: str, args: list = None):
        self.name = name
        self.command = command
        self.args = args or []
        self.finished = Signal[dict](name=f'{name}_finished')
        self._pid = None
        self._process_manager = AsyncProcessManager()

    def connectSignal(self, origin, signal_name: str, handler):
        """Connect a signal to a handler (mimics CPluginScript API)."""
        if signal_name == 'finished' and hasattr(origin, 'finished'):
            origin.finished.connect(handler, weak=False)
        else:
            raise ValueError(f"Unknown signal '{signal_name}' on {origin}")

    def process(self):
        """Start the async process."""
        # Create handler that will be called when process finishes
        handler = [self._on_process_finished, {}]

        # Start process
        self._pid = self._process_manager.startProcess(
            command=self.command,
            args=self.args,
            handler=handler,
            ifAsync=True
        )

        return MockPlugin.SUCCEEDED

    def _on_process_finished(self, pid):
        """Called when subprocess finishes."""
        # Get exit status
        exit_code = self._process_manager.getJobData(pid, 'exitCode')

        # Determine finish status
        if exit_code == 0:
            finish_status = MockPlugin.SUCCEEDED
        else:
            finish_status = MockPlugin.FAILED

        # Emit finished signal with status dict
        status_dict = {
            'pid': pid,
            'finishStatus': finish_status
        }

        self.finished.emit(status_dict)


class MockPipeline:
    """
    Mock pipeline that chains two async plugins.

    1. Run first plugin (like mtzdump)
    2. When it finishes, run second plugin (like pdbset)
    3. Report final status
    """

    SUCCEEDED = 0
    FAILED = 1

    def __init__(self, name: str):
        self.name = name
        self.finished = Signal[dict](name=f'{name}_finished')
        self.plugin1 = None
        self.plugin2 = None
        self.final_status = None

    def connectSignal(self, origin, signal_name: str, handler):
        """Connect a signal to a handler."""
        if signal_name == 'finished' and hasattr(origin, 'finished'):
            origin.finished.connect(handler, weak=False)
        else:
            raise ValueError(f"Unknown signal '{signal_name}'")

    def process(self):
        """Start the pipeline."""
        # Create first plugin (like mtzdump)
        self.plugin1 = MockPlugin(
            name='plugin1',
            command='echo',
            args=['plugin1 output']
        )

        # Connect to handler
        self.connectSignal(self.plugin1, 'finished', self.on_plugin1_finished)

        # Start it
        self.plugin1.process()

        return MockPipeline.SUCCEEDED

    @Slot(dict)
    def on_plugin1_finished(self, status_dict):
        """Called when first plugin finishes."""
        status = status_dict['finishStatus']

        if status == MockPlugin.FAILED:
            self.final_status = MockPipeline.FAILED
            self.finished.emit({'finishStatus': self.final_status})
            return

        # Create second plugin (like pdbset)
        self.plugin2 = MockPlugin(
            name='plugin2',
            command='echo',
            args=['plugin2 output']
        )

        # Connect to final handler
        self.connectSignal(self.plugin2, 'finished', self.on_plugin2_finished)

        # Start it
        self.plugin2.process()

    @Slot(dict)
    def on_plugin2_finished(self, status_dict):
        """Called when second plugin finishes."""
        status = status_dict['finishStatus']

        self.final_status = status
        self.finished.emit({'finishStatus': self.final_status})


class TestAsyncPluginIntegration:
    """Test async plugin execution with shell commands."""

    def test_single_plugin_async_execution(self):
        """Test single plugin with async subprocess."""
        # Track completion
        result = {'status': None, 'completed': False}

        @Slot(dict)
        def on_finished(status_dict):
            result['status'] = status_dict['finishStatus']
            result['completed'] = True

        # Create and run plugin
        plugin = MockPlugin(name='test_plugin', command='echo', args=['hello'])
        plugin.finished.connect(on_finished, weak=False)
        plugin.process()

        # Wait for completion
        for _ in range(50):  # Max 5 seconds
            if result['completed']:
                break
            time.sleep(0.1)

        # Verify
        assert result['completed'], "Plugin did not complete"
        assert result['status'] == MockPlugin.SUCCEEDED, "Plugin failed"

    def test_plugin_failure_handling(self):
        """Test that plugin handles subprocess failure."""
        result = {'status': None, 'completed': False}

        @Slot(dict)
        def on_finished(status_dict):
            result['status'] = status_dict['finishStatus']
            result['completed'] = True

        # Create plugin that will fail
        plugin = MockPlugin(name='fail_plugin', command='false')
        plugin.finished.connect(on_finished, weak=False)
        plugin.process()

        # Wait for completion
        for _ in range(50):
            if result['completed']:
                break
            time.sleep(0.1)

        # Verify failure was captured
        assert result['completed'], "Plugin did not complete"
        assert result['status'] == MockPlugin.FAILED, "Plugin should have failed"

    def test_pipeline_two_step_async(self):
        """Test pipeline with two async plugins in sequence."""
        result = {'status': None, 'completed': False}

        @Slot(dict)
        def on_pipeline_finished(status_dict):
            result['status'] = status_dict['finishStatus']
            result['completed'] = True

        # Create and run pipeline
        pipeline = MockPipeline(name='test_pipeline')
        pipeline.finished.connect(on_pipeline_finished, weak=False)
        pipeline.process()

        # Wait for completion (both plugins need to finish)
        for _ in range(100):  # Max 10 seconds
            if result['completed']:
                break
            time.sleep(0.1)

        # Verify both plugins ran
        assert result['completed'], "Pipeline did not complete"
        assert result['status'] == MockPipeline.SUCCEEDED, "Pipeline failed"
        assert pipeline.plugin1 is not None, "First plugin was not created"
        assert pipeline.plugin2 is not None, "Second plugin was not created"

    def test_pipeline_failure_stops_chain(self):
        """Test that pipeline stops if first plugin fails."""
        result = {'status': None, 'completed': False}

        @Slot(dict)
        def on_pipeline_finished(status_dict):
            result['status'] = status_dict['finishStatus']
            result['completed'] = True

        # Create pipeline with failing first plugin
        class FailingPipeline(MockPipeline):
            def process(self):
                # First plugin will fail
                self.plugin1 = MockPlugin(
                    name='failing_plugin',
                    command='false'  # Always fails
                )
                self.connectSignal(self.plugin1, 'finished', self.on_plugin1_finished)
                self.plugin1.process()
                return MockPipeline.SUCCEEDED

        pipeline = FailingPipeline(name='failing_pipeline')
        pipeline.finished.connect(on_pipeline_finished, weak=False)
        pipeline.process()

        # Wait for completion
        for _ in range(50):
            if result['completed']:
                break
            time.sleep(0.1)

        # Verify
        assert result['completed'], "Pipeline did not complete"
        assert result['status'] == MockPipeline.FAILED, "Pipeline should have failed"
        assert pipeline.plugin2 is None, "Second plugin should not have been created"

    def test_plugin_with_log_file(self):
        """Test plugin that writes output to log file."""
        result = {'status': None, 'completed': False}

        @Slot(dict)
        def on_finished(status_dict):
            result['status'] = status_dict['finishStatus']
            result['completed'] = True

        # Create temp directory for log
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / 'test.log'

            # Create plugin with log output
            class LoggingPlugin(MockPlugin):
                def __init__(self):
                    super().__init__(name='logging_plugin', command='echo', args=['log content'])
                    self.log_file = log_file

                def process(self):
                    handler = [self._on_process_finished, {}]
                    self._pid = self._process_manager.startProcess(
                        command=self.command,
                        args=self.args,
                        logFile=str(self.log_file),
                        handler=handler,
                        ifAsync=True
                    )
                    return MockPlugin.SUCCEEDED

            plugin = LoggingPlugin()
            plugin.finished.connect(on_finished, weak=False)
            plugin.process()

            # Wait for completion
            for _ in range(50):
                if result['completed']:
                    break
                time.sleep(0.1)

            # Verify
            assert result['completed'], "Plugin did not complete"
            assert result['status'] == MockPlugin.SUCCEEDED, "Plugin failed"
            assert log_file.exists(), "Log file was not created"

            # Check log content
            content = log_file.read_text()
            assert 'log content' in content, "Log does not contain expected output"

    def test_multiple_sequential_plugins(self):
        """Test running multiple plugins sequentially."""
        results = []
        completion_count = {'count': 0}

        def make_handler(plugin_num):
            @Slot(dict)
            def handler(status_dict):
                results.append({
                    'plugin': plugin_num,
                    'status': status_dict['finishStatus']
                })
                completion_count['count'] += 1
            return handler

        # Create and run 3 plugins in sequence
        plugins = []
        for i in range(3):
            plugin = MockPlugin(
                name=f'plugin_{i}',
                command='echo',
                args=[f'plugin {i}']
            )
            plugin.finished.connect(make_handler(i), weak=False)
            plugins.append(plugin)

        # Start them all (they'll run concurrently)
        for plugin in plugins:
            plugin.process()

        # Wait for all to complete
        for _ in range(100):
            if completion_count['count'] == 3:
                break
            time.sleep(0.1)

        # Verify all completed
        assert completion_count['count'] == 3, f"Only {completion_count['count']} plugins completed"
        assert len(results) == 3, "Not all results recorded"

        # Verify all succeeded
        for result in results:
            assert result['status'] == MockPlugin.SUCCEEDED, \
                f"Plugin {result['plugin']} failed"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
