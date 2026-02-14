"""
Test ctruncate plugin for converting intensity data to structure factors.

This test requires:
- CCP4I2_ROOT environment variable set
- CCP4 installed and ctruncate executable available
  (run via: ./run_test.sh which sets up the CCP4 environment)
"""

import pytest
import os
import shutil
from pathlib import Path
from ccp4i2.core.task_manager.plugin_registry import get_plugin_class


def get_mtz_columns(mtz_path):
    """Get column names from an MTZ file using gemmi."""
    import gemmi
    mtz = gemmi.read_mtz_file(str(mtz_path))
    return [col.label for col in mtz.columns]


def check_ccp4_available():
    """Check if CCP4 environment is available."""
    # Check if ctruncate is in PATH
    ctruncate_path = shutil.which('ctruncate')
    if ctruncate_path:
        return True

    # Check if CBIN is set
    cbin = os.environ.get('CBIN')
    if cbin and Path(cbin, 'ctruncate').exists():
        return True

    return False


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
@pytest.mark.skipif(
    not check_ccp4_available(),
    reason="CCP4 not available. Run: source /Applications/ccp4-*/bin/ccp4.setup-sh"
)
def test_ctruncate_intensity_to_fmean(tmp_path):
    """
    Test ctruncate plugin converts intensities to FMEAN (F, SIGF).

    This test uses the actual ctruncate executable to perform French-Wilson
    conversion from intensity data to structure factor amplitudes.

    Args:
        tmp_path: Pytest fixture providing a temporary directory
    """
    from ccp4i2.core.CCP4PluginScript import CPluginScript

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Get ctruncate plugin class
    ctruncate_class = get_plugin_class('ctruncate')
    assert ctruncate_class is not None, "ctruncate plugin not found"

    # Create ctruncate instance
    ctruncate = ctruncate_class(workDirectory=str(tmp_path))

    print(f"\nWork directory: {tmp_path}")
    print(f"TASKCOMMAND: {ctruncate.TASKCOMMAND}")

    # Verify ctruncate executable is accessible
    ctruncate_exe = shutil.which(ctruncate.TASKCOMMAND)
    print(f"ctruncate executable: {ctruncate_exe}")
    assert ctruncate_exe is not None, "ctruncate executable not found in PATH"

    # Set input file (intensity data)
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz"
    )
    assert Path(input_file).exists(), f"Input file not found: {input_file}"

    ctruncate.container.inputData.HKLIN.setFullPath(input_file)
    print(f"Input MTZ: {input_file}")

    # Check input columns
    input_columns = get_mtz_columns(input_file)
    print(f"Input columns: {input_columns}")

    # Configure for IMEAN -> FMEAN conversion
    # The file should have anomalous intensities (Iplus, SIGIplus, Iminus, SIGIminus)
    # IMPORTANT: Must use .set() to mark the container as set, not individual assignment
    if 'Iplus' in input_columns and 'Iminus' in input_columns:
        print("Input has anomalous intensities (Iplus/Iminus)")
        print(f"DEBUG: ISIGIanom type = {type(ctruncate.container.inputData.ISIGIanom)}")
        print(f"DEBUG: ISIGIanom class = {ctruncate.container.inputData.ISIGIanom.__class__.__name__}")
        print(f"DEBUG: ISIGIanom MRO = {[c.__name__ for c in type(ctruncate.container.inputData.ISIGIanom).__mro__[:5]]}")
        ctruncate.container.inputData.ISIGIanom.set({
            'Ip': 'Iplus',
            'SIGIp': 'SIGIplus',
            'Im': 'Iminus',
            'SIGIm': 'SIGIminus'
        })
        print(f"DEBUG: ISIGIanom.isSet() = {ctruncate.container.inputData.ISIGIanom.isSet()}")
        print(f"DEBUG: ISIGIanom.Ip type = {type(ctruncate.container.inputData.ISIGIanom.Ip)}")
        print(f"DEBUG: ISIGIanom.Ip value = {ctruncate.container.inputData.ISIGIanom.Ip}")
        print(f"DEBUG: ISIGIanom._value_states = {getattr(ctruncate.container.inputData.ISIGIanom, '_value_states', {})}")

        # Clear unused column groups to prevent ctruncate from using them
        for unused_group in ['ISIGI', 'FSIGF', 'FSIGFanom']:
            if hasattr(ctruncate.container.inputData, unused_group):
                group = getattr(ctruncate.container.inputData, unused_group)
                if hasattr(group, '_is_set'):
                    group._is_set = False
                if hasattr(group, '_column_mapping'):
                    group._column_mapping = {}

    elif 'I' in input_columns:
        print("Input has mean intensities (I/SIGI)")
        ctruncate.container.inputData.ISIGI.set({
            'I': 'I',
            'SIGI': 'SIGI'
        })
        print(f"DEBUG: ISIGI.isSet() = {ctruncate.container.inputData.ISIGI.isSet()}")

        # Clear unused column groups to prevent ctruncate from using them
        for unused_group in ['ISIGIanom', 'FSIGF', 'FSIGFanom']:
            if hasattr(ctruncate.container.inputData, unused_group):
                group = getattr(ctruncate.container.inputData, unused_group)
                if hasattr(group, '_is_set'):
                    group._is_set = False
                if hasattr(group, '_column_mapping'):
                    group._column_mapping = {}

    else:
        pytest.skip(f"Unexpected input columns: {input_columns}")

    # Set output file
    output_file = str(tmp_path / "output_fmean.mtz")
    ctruncate.container.outputData.HKLOUT.setFullPath(output_file)
    print(f"Output MTZ: {output_file}")

    # Configure to output mini-MTZ with FMEAN format
    ctruncate.container.controlParameters.OUTPUTMINIMTZ.set(True)
    ctruncate.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG.set(4)  # FMEAN

    # Run ctruncate (this will call makeCommandAndScript internally)
    print("\n" + "="*60)
    print("Running ctruncate...")
    print("="*60)

    status = ctruncate.process()

    # Show command line that was generated
    print(f"\nCommand line that was executed:")
    print(f"  {ctruncate.TASKCOMMAND} {' '.join(ctruncate.commandLine)}")

    print(f"\nProcess status: {status}")
    print(f"SUCCEEDED = {CPluginScript.SUCCEEDED}")
    print(f"RUNNING = {CPluginScript.RUNNING}")
    print(f"FAILED = {CPluginScript.FAILED}")

    # Check for errors
    error_count = ctruncate.errorReport.count()
    if error_count > 0:
        print(f"\nErrors ({error_count}):")
        print(ctruncate.errorReport.report())

    # List files created
    created_files = list(tmp_path.iterdir())
    print(f"\nFiles created ({len(created_files)}):")
    for f in created_files:
        if f.is_file():
            size_kb = f.stat().st_size / 1024
            print(f"  - {f.name} ({size_kb:.1f} KB)")

    # Check if output file was created
    assert Path(output_file).exists(), f"Output file not created: {output_file}"
    assert Path(output_file).stat().st_size > 0, "Output file is empty"

    print(f"\n✅ Output file created: {output_file}")

    # Verify output has structure factor columns
    output_columns = get_mtz_columns(output_file)
    print(f"Output columns: {output_columns}")

    # ctruncate outputs FMEAN/SIGFMEAN (or sometimes F/SIGF depending on mode)
    has_f = 'F' in output_columns or 'FMEAN' in output_columns
    has_sigf = 'SIGF' in output_columns or 'SIGFMEAN' in output_columns

    assert has_f, f"F or FMEAN column not found in output. Columns: {output_columns}"
    assert has_sigf, f"SIGF or SIGFMEAN column not found in output. Columns: {output_columns}"

    print("✅ Output contains structure factor columns (French-Wilson conversion successful)")

    # Verify process completed successfully
    assert status == CPluginScript.SUCCEEDED or status == CPluginScript.RUNNING, \
        f"Process did not complete successfully. Status: {status}"


if __name__ == "__main__":
    # Allow running this test directly
    import sys
    sys.exit(pytest.main([__file__, "-v", "-s"]))
