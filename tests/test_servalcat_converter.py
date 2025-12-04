"""
Unit tests for servalcat fw French-Wilson converter.

Tests the new ServalcatConverter that replaces ctruncate for
intensity-to-amplitude conversions.
"""
import pytest
from pathlib import Path


def test_servalcat_fw_available():
    """Check that servalcat fw is available in the environment."""
    import subprocess
    result = subprocess.run(
        ['which', 'servalcat'], capture_output=True, text=True)
    assert result.returncode == 0, (
        "servalcat not found in PATH. Source CCP4 environment.")


def test_imean_to_fmean_conversion(tmp_path):
    """Test IMEAN → FMEAN conversion using servalcat fw."""
    from core.CCP4XtalData import CObsDataFile
    import shutil

    # Use existing test data (IMEAN from gamma dataset)
    test_data = (
        Path(__file__).parent.parent / 'demo_data' / 'gamma' /
        'HKLOUT_unmerged.mtz')

    if not test_data.exists():
        pytest.skip(f"Test data not found: {test_data}")

    # Create a temporary copy to work with
    input_mtz = tmp_path / "input_intensities.mtz"
    shutil.copy2(test_data, input_mtz)

    # Create CObsDataFile instance
    obs_file = CObsDataFile()
    obs_file.setFullPath(str(input_mtz))

    # Auto-detect content flag
    obs_file.setContentFlag()

    # Test expects IMEAN data
    if int(obs_file.contentFlag) != obs_file.CONTENT_FLAG_IMEAN:
        pytest.skip(
            f"Test data is not IMEAN (contentFlag="
            f"{obs_file.contentFlag})")

    # Test conversion
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    output_path = obs_file.as_FMEAN(work_directory=str(work_dir))

    # Verify output exists
    assert Path(output_path).exists(), (
        f"Output file not created: {output_path}")

    # Verify it's an MTZ file
    import gemmi
    mtz = gemmi.read_mtz_file(output_path)

    # Check for FMEAN columns
    column_labels = [col.label for col in mtz.columns]
    assert 'F' in column_labels, "F column not found in output"
    assert 'SIGF' in column_labels, "SIGF column not found in output"

    # Check servalcat subdirectory was created
    servalcat_dir = work_dir / 'servalcat'
    assert servalcat_dir.exists(), "servalcat subdirectory not created"

    # Check monolithic MTZ exists (for caching)
    monolithic_mtz_pattern = list(servalcat_dir.glob("*.mtz"))
    assert len(monolithic_mtz_pattern) > 0, (
        "Monolithic MTZ not found in servalcat directory")

    print(f"✅ IMEAN → FMEAN conversion successful: {output_path}")
    print(f"   Columns: {column_labels}")
    print(f"   Monolithic MTZ cached: {monolithic_mtz_pattern[0]}")


def test_monolithic_mtz_contains_all_columns(tmp_path):
    """Test that servalcat creates monolithic MTZ with all outputs."""
    from core.CCP4XtalData import CObsDataFile
    import shutil

    # Use IMEAN test data from gamma dataset
    test_data = (
        Path(__file__).parent.parent / 'demo_data' / 'gamma' /
        'HKLOUT_unmerged.mtz')

    if not test_data.exists():
        pytest.skip(f"Test data not found: {test_data}")

    input_mtz = tmp_path / "input_intensities.mtz"
    shutil.copy2(test_data, input_mtz)

    obs_file = CObsDataFile()
    obs_file.setFullPath(str(input_mtz))
    obs_file.setContentFlag()

    if int(obs_file.contentFlag) != obs_file.CONTENT_FLAG_IMEAN:
        pytest.skip("Test data is not IMEAN")

    work_dir = tmp_path / "work"
    work_dir.mkdir()

    # Run conversion (IMEAN → FMEAN)
    obs_file.as_FMEAN(work_directory=str(work_dir))

    # Check monolithic MTZ
    servalcat_dir = work_dir / 'servalcat'
    monolithic_mtz = list(servalcat_dir.glob("*.mtz"))[0]

    # Verify monolithic MTZ has expected columns
    import gemmi
    mtz = gemmi.read_mtz_file(str(monolithic_mtz))
    column_labels = [col.label for col in mtz.columns]

    # For IMEAN input, servalcat fw produces:
    # F, SIGF (FMEAN)
    # I, SIGI (IMEAN)

    print(f"✅ Monolithic MTZ columns: {column_labels}")

    # Should contain FMEAN columns
    assert 'F' in column_labels, "F not in monolithic MTZ"
    assert 'SIGF' in column_labels, "SIGF not in monolithic MTZ"

    # Should contain IMEAN columns
    assert 'I' in column_labels, "I not in monolithic MTZ"
    assert 'SIGI' in column_labels, "SIGI not in monolithic MTZ"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
