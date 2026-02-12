#!/usr/bin/env python3
"""Test freerflag saveParams() behavior."""

import tempfile
from pathlib import Path

from ccp4i2.core.task_manager.plugin_registry import get_plugin_class

def test_freerflag_save():
    """Create a freerflag plugin and save params to see what gets written."""

    # Create a temporary directory for the test
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)

        # Get the plugin class
        plugin_class = get_plugin_class('freerflag')

        # Create plugin instance
        print(f"\n=== Creating plugin ===")
        plugin = plugin_class(parent=None, name='freerflag_test')

        # Check RESMAX state immediately after plugin creation
        print(f"\n=== Immediately After Plugin Creation ===")
        resmax_fresh = plugin.container.controlParameters.RESMAX
        print(f"RESMAX id: {id(resmax_fresh)}")
        print(f"RESMAX value_state: {resmax_fresh.getValueState('value')}")
        print(f"RESMAX _value_states: {resmax_fresh._value_states}")

        # Check which RESMAX object IDs were created during def_xml parsing
        print(f"\nLooking back at debug output for all RESMAX object IDs...")

        # Set the work directory
        plugin.workDirectory = workdir

        # Check again after setting workDirectory
        print(f"\n=== After Setting workDirectory ===")
        print(f"RESMAX value_state: {plugin.container.controlParameters.RESMAX.getValueState('value')}")

        print(f"\n=== Initial State ===")
        resmax = plugin.container.controlParameters.RESMAX
        print(f"RESMAX value: {resmax.value}")
        print(f"RESMAX value_state: {resmax.getValueState('value')}")
        print(f"RESMAX isSet(): {resmax.isSet()}")
        print(f"RESMAX isSet(allowDefault=False): {resmax.isSet(allowDefault=False)}")
        print(f"RESMAX object id: {id(resmax)}")
        print(f"RESMAX has default qualifier: {resmax.get_qualifier('default')}")

        # Check _value_states directly
        if hasattr(resmax, '_value_states'):
            print(f"RESMAX _value_states dict: {resmax._value_states}")

        # Try to save params
        print(f"\n=== Right Before saveParams() ===")
        resmax_before = plugin.container.controlParameters.RESMAX
        print(f"RESMAX value_state: {resmax_before.getValueState('value')}")
        print(f"RESMAX object id: {id(resmax_before)}")
        print(f"Same object? {resmax is resmax_before}")

        print(f"\n=== Calling saveParams() ===")
        output_file = workdir / "test_params.xml"
        error = plugin.saveDataToXml(str(output_file))

        if error:
            print(f"Error saving: {error}")
        else:
            print(f"Saved successfully to: {output_file}")

        # Read and display the saved XML
        if output_file.exists():
            print(f"\n=== Saved XML Content ===")
            with open(output_file, 'r') as f:
                content = f.read()
                print(content)

            # Check if RESMAX appears in the XML (look for the exact tag)
            if '<RESMAX>' in content:
                print("\n❌ RESMAX appears in saved XML (should be excluded!)")
            else:
                print("\n✅ RESMAX correctly excluded from saved XML")
        else:
            print(f"\n❌ Output file not created: {output_file}")

if __name__ == '__main__':
    test_freerflag_save()
