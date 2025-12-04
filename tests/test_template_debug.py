"""Debug template expansion for mtzdump."""

import os
import sys
import tempfile
from pathlib import Path

# Set up environment
os.environ['CCP4I2_ROOT'] = '/Users/nmemn/Developer/ccp4i2'
sys.path.insert(0, '/Users/nmemn/Developer/ccp4i2')

# Import mt zdump plugin
from wrappers.mtzdump.script.mtzdump import mtzdump

def test_template_expansion():
    """Test if template expansion works for mtzdump."""

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a dummy MTZ file path (doesn't need to exist for this test)
        test_mtz = Path(tmpdir) / "test.mtz"
        test_mtz.touch()  # Create empty file

        print("\n" + "="*70)
        print("Testing mtzdump template expansion")
        print("="*70)

        # Create mtzdump instance
        plugin = mtzdump(workDirectory=tmpdir, name="test_mtzdump")

        print(f"\nCOMLINETEMPLATE: {plugin.COMLINETEMPLATE}")
        print(f"COMTEMPLATE: {plugin.COMTEMPLATE}")

        # Set input file
        plugin.container.inputData.HKLIN.setFullPath(str(test_mtz))
        print(f"\nHKLIN path set to: {plugin.container.inputData.HKLIN.getFullPath()}")

        # Test find() method
        print(f"\nTesting find() method:")
        found = plugin.container.find("HKLIN")
        print(f"  container.find('HKLIN'): {found}")
        if found:
            print(f"  Type: {type(found).__name__}")
            print(f"  str(found): {str(found)}")

        found2 = plugin.container.find("inputData.HKLIN")
        print(f"  container.find('inputData.HKLIN'): {found2}")
        if found2:
            print(f"  Type: {type(found2).__name__}")
            print(f"  str(found2): {str(found2)}")

        # Try makeCommandAndScript
        print(f"\nCalling makeCommandAndScript()...")
        error = plugin.makeCommandAndScript()
        print(f"Errors: {error}")
        print(f"Command line: {plugin.commandLine}")
        print(f"Command script: {plugin.commandScript}")

        print("\n" + "="*70)

if __name__ == '__main__':
    test_template_expansion()
