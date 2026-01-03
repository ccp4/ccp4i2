#!/usr/bin/env python3
"""Simple test to verify RESMAX stays NOT_SET and is excluded from XML."""

import sys
import tempfile
from pathlib import Path

from ccp4i2.core.CCP4TaskManager import TASKMANAGER

def test_resmax_not_set():
    """Verify RESMAX without default stays NOT_SET and is excluded from XML."""
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)

        # Create freerflag plugin
        plugin_class = TASKMANAGER().get_plugin_class('freerflag')
        plugin = plugin_class(parent=None, name='freerflag_test')
        plugin.workDirectory = workdir

        # Check RESMAX state
        resmax = plugin.container.controlParameters.RESMAX
        state = resmax.getValueState('value')
        print(f"RESMAX value_state: {state}")
        print(f"RESMAX value: {resmax.value}")
        print(f"RESMAX isSet(allowDefault=False): {resmax.isSet(allowDefault=False)}")

        # Save to XML
        output_file = workdir / "test_params.xml"
        plugin.saveDataToXml(str(output_file))

        # Check if RESMAX is in the XML
        with open(output_file, 'r') as f:
            content = f.read()

        if '<RESMAX>' in content:
            print("\n❌ FAIL: RESMAX appears in XML (should be excluded)")
            print(content)
            return False
        else:
            print("\n✅ PASS: RESMAX correctly excluded from XML")
            return True

if __name__ == '__main__':
    success = test_resmax_not_set()
    sys.exit(0 if success else 1)
