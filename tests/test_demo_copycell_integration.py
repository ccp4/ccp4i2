"""
Integration test for real demo_copycell plugin using CCP4 demo data.

This tests the complete async infrastructure with a real CCP4i2 plugin
and actual crystallographic data files.
"""

import pytest
import time
import tempfile
import os
import sys
from pathlib import Path
import shutil

# Ensure CCP4I2_ROOT is set
CCP4I2_ROOT = os.environ.get('CCP4I2_ROOT')
if not CCP4I2_ROOT:
    pytest.skip("CCP4I2_ROOT not set", allow_module_level=True)

# Check if demo_data exists
demo_data_path = Path(CCP4I2_ROOT) / 'demo_data' / 'mdm2'
if not demo_data_path.exists():
    pytest.skip(f"Demo data not found: {demo_data_path}", allow_module_level=True)

# Add legacy ccp4i2 to path (at the end to avoid conflicts)
if CCP4I2_ROOT not in sys.path:
    sys.path.append(CCP4I2_ROOT)

# Now try to import the real demo_copycell plugin
try:
    # Import our modern CPluginScript first
    from ccp4i2.core.CCP4PluginScript import CPluginScript as ModernCPluginScript
    from ccp4i2.core.base_object.signal_system import Slot

    # Try to import the legacy plugin
    # This might fail if Qt dependencies aren't satisfied
    from ccp4i2.wrappers2.demo_copycell.script.demo_copycell import demo_copycell
    DEMO_COPYCELL_AVAILABLE = True
except ImportError as e:
    print(f"Could not import demo_copycell: {e}")
    DEMO_COPYCELL_AVAILABLE = False


class TestDemoCopycellIntegration:
    """Integration test for real demo_copycell plugin."""

    @pytest.mark.skipif(not DEMO_COPYCELL_AVAILABLE,
                        reason="demo_copycell plugin not available")
    def test_demo_copycell_with_mdm2_data(self):
        """
        Run real demo_copycell plugin with MDM2 test data.

        This test:
        1. Copies cell parameters from mdm2_unmerged.mtz
        2. Applies them to 4qo4.pdb using real mtzdump and pdbset
        3. Verifies async execution completes successfully
        4. Checks output PDB file is created

        Expected behavior:
        - mtzdump extracts cell: ~42.7 42.7 42.7 90 90 90 (P213)
        - pdbset writes output PDB with updated cell
        """
        # Input files from demo_data
        input_mtz = demo_data_path / 'mdm2_unmerged.mtz'
        input_pdb = demo_data_path / '4qo4.pdb'

        if not input_mtz.exists():
            pytest.skip(f"Test data not found: {input_mtz}")
        if not input_pdb.exists():
            pytest.skip(f"Test data not found: {input_pdb}")

        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy input files to work directory
            work_mtz = Path(tmpdir) / 'input.mtz'
            work_pdb = Path(tmpdir) / 'input.pdb'
            output_pdb = Path(tmpdir) / 'output.pdb'

            shutil.copy(input_mtz, work_mtz)
            shutil.copy(input_pdb, work_pdb)

            print(f"\n{'='*70}")
            print("Running demo_copycell with MDM2 data")
            print(f"{'='*70}")
            print(f"Input MTZ: {work_mtz}")
            print(f"Input PDB: {work_pdb}")
            print(f"Output PDB: {output_pdb}")
            print(f"{'='*70}\n")

            # Track completion
            result = {
                'status': None,
                'completed': False,
                'mtzdump_cell': None,
                'error': None
            }

            @Slot(dict)
            def on_finished(status_dict):
                """Handler for pipeline completion."""
                result['status'] = status_dict.get('finishStatus')
                result['completed'] = True
                print(f"\n{'='*70}")
                print("Pipeline finished!")
                print(f"Status: {result['status']}")
                print(f"{'='*70}\n")

            try:
                # Create demo_copycell instance
                pipeline = demo_copycell(
                    workDirectory=tmpdir,
                    name="test_demo_copycell_mdm2"
                )

                # Set input files using setFullPath API
                pipeline.container.inputData.HKLIN.setFullPath(str(work_mtz))
                pipeline.container.inputData.XYZIN.setFullPath(str(work_pdb))
                pipeline.container.outputData.XYZOUT.setFullPath(str(output_pdb))

                # Connect finished signal
                pipeline.finished.connect(on_finished, weak=False)

                # Start pipeline
                print("Starting demo_copycell pipeline...")
                pipeline.process()

                # Wait for completion
                print("Waiting for completion...\n")
                max_wait = 600  # 60 seconds max (real programs take time)
                for i in range(max_wait):
                    if result['completed']:
                        break
                    time.sleep(0.1)
                    if (i + 1) % 50 == 0:
                        print(f"  ... waiting ({(i+1)/10:.1f}s)")

                # Check if we timed out
                if not result['completed']:
                    pytest.fail("Pipeline did not complete within 60 seconds")

                # Verify completion
                print(f"\n{'='*70}")
                print("Verification Results")
                print(f"{'='*70}")
                print(f"Completed: {result['completed']}")
                print(f"Status: {result['status']}")

                # Check status
                assert result['status'] == ModernCPluginScript.SUCCEEDED, \
                    f"Pipeline failed with status {result['status']}"

                # Verify sub-plugins were created
                assert hasattr(pipeline, 'mtzdump'), \
                    "mtzdump plugin was not created"
                assert hasattr(pipeline, 'pdbset'), \
                    "pdbset plugin was not created"

                print(f"mtzdump plugin: {pipeline.mtzdump}")
                print(f"pdbset plugin: {pipeline.pdbset}")

                # Try to get cell from mtzdump output
                if hasattr(pipeline.mtzdump, 'container'):
                    if hasattr(pipeline.mtzdump.container, 'outputData'):
                        if hasattr(pipeline.mtzdump.container.outputData, 'CELL'):
                            cell = pipeline.mtzdump.container.outputData.CELL
                            print(f"Cell extracted by mtzdump: {cell}")
                            result['mtzdump_cell'] = cell

                # Verify output file
                print(f"Output PDB exists: {output_pdb.exists()}")
                if output_pdb.exists():
                    print(f"Output PDB size: {output_pdb.stat().st_size} bytes")

                assert output_pdb.exists(), \
                    "Output PDB file was not created"
                assert output_pdb.stat().st_size > 0, \
                    "Output PDB file is empty"

                # Read first few lines of output to verify CRYST1 record
                with open(output_pdb, 'r') as f:
                    for line in f:
                        if line.startswith('CRYST1'):
                            print(f"CRYST1 record: {line.strip()}")
                            break

                print(f"\n✅ All checks passed!")
                print(f"{'='*70}\n")

            except Exception as e:
                result['error'] = str(e)
                print(f"\n❌ Error during pipeline execution: {e}")
                import traceback
                traceback.print_exc()
                pytest.fail(f"Pipeline execution failed: {e}")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
