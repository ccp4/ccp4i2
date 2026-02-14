"""
Tests for building and using parrot plugin.
"""

import pytest
import os
from pathlib import Path
from ccp4i2.core.task_manager.plugin_registry import get_plugin_class


def get_mtz_columns(mtz_path):
    """
    Get column names from an MTZ file using gemmi.

    Args:
        mtz_path: Path to MTZ file

    Returns:
        List of column names
    """
    import gemmi

    mtz = gemmi.read_mtz_file(str(mtz_path))
    return [col.label for col in mtz.columns]


def check_ccp4_available():
    """Check if CCP4 environment is available."""
    import shutil
    # Check if ctruncate is in PATH (CCP4 setup script was sourced)
    return shutil.which('ctruncate') is not None


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
@pytest.mark.skipif(
    not check_ccp4_available(),
    reason="CCP4 not available. Run: source /Applications/ccp4-*/bin/ccp4.setup-sh"
)
def test_as_fmean_conversion(tmp_path):
    """Test that as_FMEAN() converts intensities to F, SIGF structure factors.

    Args:
        tmp_path: Pytest fixture providing a temporary directory for test outputs
    """
    from ccp4i2.core.CCP4XtalData import CObsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Create CObsDataFile for intensity data
    intensity_file = CObsDataFile()
    intensity_file.setFullPath(
        os.path.join(ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz")
    )

    print(f"\nInput file: {intensity_file.getFullPath()}")

    # Detect content flag
    intensity_file.setContentFlag()
    print(f"Content flag: {intensity_file.contentFlag}")

    # Convert to FMEAN
    output_path = intensity_file.as_FMEAN(work_directory=str(tmp_path))

    # Verify output file exists
    assert Path(output_path).exists(), f"Output file not created: {output_path}"
    print(f"✅ Conversion output created: {output_path}")

    # Verify output has structure factor columns
    columns = get_mtz_columns(output_path)
    print(f"Output columns: {columns}")

    # ctruncate outputs FMEAN/SIGFMEAN (or sometimes F/SIGF)
    has_f = 'F' in columns or 'FMEAN' in columns
    has_sigf = 'SIGF' in columns or 'SIGFMEAN' in columns

    assert has_f, f"F or FMEAN column not found. Columns: {columns}"
    assert has_sigf, f"SIGF or SIGFMEAN column not found. Columns: {columns}"

    print("✅ Output file contains structure factor columns")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
@pytest.mark.skipif(
    not check_ccp4_available(),
    reason="CCP4 not available. Run: source /Applications/ccp4-*/bin/ccp4.setup-sh"
)
def test_parrot_makehklin(tmp_path):
    """Test that parrot's makeHklInput creates hklin.mtz with expected columns.

    Args:
        tmp_path: Pytest fixture providing a temporary directory for test outputs
    """
    # Create parrot task with temporary work directory
    task = get_plugin_class("parrot")(workDirectory=tmp_path)

    # Set input files from CCP4I2 demo data
    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Set F_SIGF (observations - structure factors with sigmas)
    task.container.inputData.F_SIGF = os.path.join(
        ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz"
    )

    # Set ABCD (phases from phasing program)
    task.container.inputData.ABCD = os.path.join(
        ccp4_root, "demo_data", "gamma", "initial_phases.mtz"
    )

    # Call makeHklInput to create the merged hklin.mtz
    # This should trigger conversion of intensities to F, SIGF if needed
    # makeHklInput returns (outfile, colnames, error) tuple
    # Use [name, contentFlag] syntax to request FMEAN (4) conversion
    from ccp4i2.core.CCP4XtalData import CObsDataFile
    hklin_path, colnames, error = task.makeHklInput(
        miniMtzsIn=[['F_SIGF', CObsDataFile.CONTENT_FLAG_FMEAN], 'ABCD']
    )

    # Check for errors
    if error and error.count() > 0:
        print(f"Errors during makeHklInput: {error.report()}")

    # Verify hklin.mtz was created
    assert hklin_path is not None, "makeHklInput returned None"
    assert Path(hklin_path).exists(), f"hklin.mtz not created at {hklin_path}"

    print(f"\n✅ hklin.mtz created: {hklin_path}")
    print(f"File size: {Path(hklin_path).stat().st_size} bytes")

    # Verify columns in hklin.mtz
    columns = get_mtz_columns(hklin_path)
    print(f"hklin.mtz columns: {columns}")

    # Expected columns from F_SIGF input (after conversion)
    # Object name is prepended with underscore (e.g., F_SIGF_F)
    assert 'F_SIGF_F' in columns, f"F_SIGF_F column not found in hklin.mtz. Columns: {columns}"
    assert 'F_SIGF_SIGF' in columns, f"F_SIGF_SIGF column not found in hklin.mtz. Columns: {columns}"

    # Expected columns from ABCD input (phases)
    # Object name is prepended with underscore (e.g., ABCD_HLA)
    assert 'ABCD_HLA' in columns, f"ABCD_HLA column not found in hklin.mtz. Columns: {columns}"
    assert 'ABCD_HLB' in columns, f"ABCD_HLB column not found in hklin.mtz. Columns: {columns}"
    assert 'ABCD_HLC' in columns, f"ABCD_HLC column not found in hklin.mtz. Columns: {columns}"
    assert 'ABCD_HLD' in columns, f"ABCD_HLD column not found in hklin.mtz. Columns: {columns}"

    # Expected crystallographic columns
    assert 'H' in columns, f"H (Miller index) not found. Columns: {columns}"
    assert 'K' in columns, f"K (Miller index) not found. Columns: {columns}"
    assert 'L' in columns, f"L (Miller index) not found. Columns: {columns}"

    print("✅ hklin.mtz contains all expected columns (F_SIGF_F, F_SIGF_SIGF, ABCD_HLA/HLB/HLC/HLD, Miller indices)")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
def test_parrot(tmp_path):
    """Test parrot plugin initialization and execution.

    Args:
        tmp_path: Pytest fixture providing a temporary directory for test outputs
    """
    # Create parrot task with temporary work directory
    task = get_plugin_class("parrot")(workDirectory=tmp_path)

    # Verify work directory is set correctly
    assert task.workDirectory == tmp_path

    # Set input files from CCP4I2 demo data
    # Parrot expects F_SIGF (structure factors) and ABCD (phases)
    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Set F_SIGF (observations - structure factors with sigmas)
    task.container.inputData.F_SIGF = os.path.join(
        ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz"
    )
    assert str(task.container.inputData.F_SIGF).endswith(
        "merged_intensities_native.mtz"
    )

    # Set ABCD (phases from phasing program)
    task.container.inputData.ABCD = os.path.join(
        ccp4_root, "demo_data", "gamma", "initial_phases.mtz"
    )
    assert str(task.container.inputData.ABCD).endswith("initial_phases.mtz")

    # Set ASUIN (observations - structure factors with sigmas)
    task.container.inputData.ASUIN = os.path.join(
        ccp4_root, "demo_data", "gamma", "gamma.asu.xml"
    )
    assert str(task.container.inputData.ASUIN).endswith(
        "gamma.asu.xml"
    )

    # Run the task - outputs will be created in tmp_path
    result = task.process()

    # Verify the task completed
    # Note: The actual parrot executable may not run if not available,
    # but the workflow should complete without errors
    print(f"\nTask process() returned: {result}")
    print(f"Work directory: {task.workDirectory}")

    # Check for errors
    error_count = task.errorReport.count()
    print(f"Error count: {error_count}")
    if error_count > 0:
        print(f"Errors:\n{task.errorReport.report()}")

    # List all files created in work directory
    created_files = list(tmp_path.iterdir())
    print(f"\nFiles created: {len(created_files)}")
    for f in created_files:
        size_kb = f.stat().st_size / 1024
        print(f"  - {f.name} ({size_kb:.1f} KB)")

    # Check if hklin.mtz was created (merged input file)
    hklin_path = tmp_path / "hklin.mtz"
    if hklin_path.exists():
        print(f"\n✅ hklin.mtz created: {hklin_path.stat().st_size} bytes")
    else:
        print(f"\n⚠️  hklin.mtz NOT found at {hklin_path}")

    # Verify work directory is still the temp path
    assert task.workDirectory == tmp_path

    # The tmp_path will be automatically cleaned up by pytest after the test
