"""
Test demo_copycell pattern with legacy int signature support.

This tests that legacy plugins with @QtCore.Slot(int) decorators work
with our new async infrastructure by using signature adaptation.
"""

import pytest
import time
import tempfile
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.signal_system import Slot


class LegacyAsyncPlugin(CPluginScript):
    """
    Plugin that mimics legacy behavior with int signature.

    This tests the signature adaptation in connectSignal().
    """

    TASKMODULE = 'test'
    TASKTITLE = 'Legacy Async Plugin'
    TASKNAME = 'legacy_async'
    TASKCOMMAND = 'echo'
    TASKVERSION = '1.0'
    ASYNCHRONOUS = True

    def __init__(self, message: str = "legacy output", **kwargs):
        super().__init__(**kwargs)
        self.message = message

    def process(self):
        """Start the async process."""
        self.commandLine = [self.message]
        return super().process()


class LegacyPipeline(CPluginScript):
    """
    Pipeline that mimics demo_copycell with legacy int signatures.

    The handlers use type hints (int) to test signature adaptation.
    """

    TASKMODULE = 'test'
    TASKTITLE = 'Legacy Pipeline Test'
    TASKNAME = 'legacy_pipeline'
    TASKVERSION = '1.0'
    ASYNCHRONOUS = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def process(self):
        """Start the pipeline."""
        print("\n" + "="*60)
        print("Starting legacy pipeline (int signatures)")
        print("="*60)

        # Create first plugin
        self.plugin1 = LegacyAsyncPlugin(
            message="plugin1 output",
            name="plugin1",
            workDirectory=self.workDirectory
        )

        # Connect with legacy int signature - note the type hint!
        self.connectSignal(self.plugin1, 'finished', self.on_plugin1_finished)

        # Start plugin1
        print("Starting plugin1...")
        self.plugin1.process()

        return self.SUCCEEDED

    def on_plugin1_finished(self, status: int):
        """
        Legacy handler with INT signature (not dict).

        This tests that connectSignal() adapts the dict emission
        to just pass the int status.
        """
        print(f"\n✅ plugin1 finished with status: {status} (type: {type(status).__name__})")
        print(f"   Status is int: {isinstance(status, int)}")

        # Verify we received an int, not a dict
        assert isinstance(status, int), f"Expected int, got {type(status).__name__}"

        if status == self.FAILED:
            self.reportStatus(self.FAILED)
            return

        # Create second plugin
        self.plugin2 = LegacyAsyncPlugin(
            message="plugin2 output",
            name="plugin2",
            workDirectory=self.workDirectory
        )

        # Connect to final handler (also with int signature)
        self.connectSignal(self.plugin2, 'finished', self.on_pipeline_finished)

        # Start plugin2
        print("Starting plugin2...")
        self.plugin2.process()

    def on_pipeline_finished(self, status: int):
        """Final handler with int signature."""
        print(f"\n✅ plugin2 finished with status: {status} (type: {type(status).__name__})")
        print(f"   Status is int: {isinstance(status, int)}")

        assert isinstance(status, int), f"Expected int, got {type(status).__name__}"
        self.reportStatus(status)


class TestDemoCopycellReal:
    """Test demo_copycell pattern with legacy int signature."""

    def test_legacy_int_signature_pipeline(self):
        """
        Test pipeline with legacy int signature handlers.

        This verifies that our connectSignal() correctly adapts:
        - Modern emission: finished.emit({'finishStatus': status})
        - Legacy handler: def on_finished(self, status: int)

        The adaptation happens via type hint inspection.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            result = {'status': None, 'completed': False}

            @Slot(dict)
            def on_finished(status_dict):
                """Modern handler that receives dict."""
                print(f"\n{'='*60}")
                print("Pipeline completed!")
                print(f"{'='*60}")
                print(f"Received dict: {status_dict}")
                result['status'] = status_dict.get('finishStatus')
                result['completed'] = True

            # Create legacy pipeline
            pipeline = LegacyPipeline(
                workDirectory=tmpdir,
                name="test_legacy_pipeline"
            )

            # Connect with modern dict handler
            pipeline.finished.connect(on_finished, weak=False)

            # Start pipeline
            status = pipeline.process()
            assert status == CPluginScript.SUCCEEDED

            # Wait for completion
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

            print(f"\n✅ All checks passed! Legacy int signature support works!")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
