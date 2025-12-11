#!/usr/bin/env python3
"""
Quick test to verify that CPluginScript.process() forwards kwargs to startProcess()
"""

import sys
import os
from pathlib import Path

# Set up environment using relative path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
os.environ['CCP4I2_ROOT'] = str(PROJECT_ROOT)
os.environ['DJANGO_SETTINGS_MODULE'] = 'ccp4x.settings'

# Add project root to path
sys.path.insert(0, str(PROJECT_ROOT))

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Modules


class TestPlugin(CPluginScript):
    """Test plugin to verify kwargs forwarding"""

    TASKNAME = 'test_kwargs_plugin'
    TASKMODULE = 'test'
    TASKTITLE = 'Test Kwargs Plugin'

    def __init__(self):
        super().__init__()
        self.received_kwargs = {}

    def checkInputData(self):
        """Skip input validation for test"""
        return None

    def checkOutputData(self):
        """Skip output setup for test"""
        return None

    def processInputFiles(self):
        """Skip input processing for test"""
        return self.SUCCEEDED

    def makeCommandAndScript(self):
        """Skip command generation for test"""
        return None

    def startProcess(self, dummy=None, **kw):
        """
        Capture kwargs that were forwarded from process()
        """
        print(f"startProcess called with kwargs: {kw}")
        self.received_kwargs = kw
        return self.SUCCEEDED


def test_kwargs_forwarding():
    """Test that process() forwards kwargs to startProcess()"""

    plugin = TestPlugin()

    # Call process() with custom kwargs
    result = plugin.process(filename="/path/to/file.mtz", pxdname="MyDataset", threshold=0.8)

    # Verify kwargs were forwarded
    assert result == plugin.SUCCEEDED, "Process should succeed"
    assert 'filename' in plugin.received_kwargs, "filename should be forwarded"
    assert 'pxdname' in plugin.received_kwargs, "pxdname should be forwarded"
    assert 'threshold' in plugin.received_kwargs, "threshold should be forwarded"
    assert plugin.received_kwargs['filename'] == "/path/to/file.mtz", "filename value should match"
    assert plugin.received_kwargs['pxdname'] == "MyDataset", "pxdname value should match"
    assert plugin.received_kwargs['threshold'] == 0.8, "threshold value should match"

    print("\n✅ SUCCESS: All kwargs were correctly forwarded from process() to startProcess()")
    print(f"   Received kwargs: {plugin.received_kwargs}")
    return True


if __name__ == "__main__":
    try:
        test_kwargs_forwarding()
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
