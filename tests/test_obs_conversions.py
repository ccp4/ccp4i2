"""
Test CObsDataFile contentFlag conversions.

This test suite verifies all supported conversions between
MTZ observation data formats:
- IPAIR (1): Anomalous intensities (I+, SIGI+, I-, SIGI-)
- FPAIR (2): Anomalous structure factors (F+, SIGF+, F-, SIGF-)
- IMEAN (3): Mean intensities (I, SIGI)
- FMEAN (4): Mean structure factors (F, SIGF)

Conversion matrix (from original CCP4XtalData.py):
                    TO
          IPAIR FPAIR IMEAN FMEAN
FROM IPAIR   0    2    2     2     (0=same, 2=ctruncate, 4=not possible)
     FPAIR   4    0    4     2
     IMEAN   4    4    0     2
     FMEAN   4    4    4     0

Requires:
- CCP4I2_ROOT environment variable
- CCP4 installed with ctruncate executable
"""

import pytest
import os
import shutil
from pathlib import Path
from ccp4i2.core.CCP4TaskManager import TASKMANAGER


def get_mtz_columns(mtz_path):
    """Get column names from an MTZ file using gemmi."""
    import gemmi
    mtz = gemmi.read_mtz_file(str(mtz_path))
    return [col.label for col in mtz.columns]


def check_ccp4_available():
    """Check if CCP4 environment is available."""
    ctruncate_path = shutil.which('ctruncate')
    if ctruncate_path:
        return True

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
def test_ipair_to_fpair_conversion(tmp_path):
    """
    Test IPAIR → FPAIR conversion.

    Converts anomalous intensities (I+/I-) to anomalous structure factor
    amplitudes (F+/F-) using French-Wilson conversion.
    """
    from ccp4i2.core.CCP4XtalData import CObsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Load IPAIR file (anomalous intensities)
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz"
    )
    assert Path(input_file).exists(), f"Input file not found: {input_file}"

    # Create file object
    ipair_file = CObsDataFile()
    ipair_file.setFullPath(input_file)
    ipair_file.setContentFlag()

    print(f"\nInput file: {input_file}")
    print(f"Input contentFlag: {ipair_file.contentFlag}")
    print(f"Input columns: {get_mtz_columns(input_file)}")

    assert int(ipair_file.contentFlag) == CObsDataFile.CONTENT_FLAG_IPAIR, \
        f"Expected IPAIR (1), got {ipair_file.contentFlag}"

    # Convert to FPAIR
    print("\n" + "=" * 60)
    print("Converting IPAIR → FPAIR...")
    print("=" * 60)

    output_file = ipair_file.as_FPAIR(work_directory=str(tmp_path))

    print(f"\n✅ Conversion completed: {output_file}")

    # Verify output exists and has correct columns
    assert Path(output_file).exists(), f"Output file not created: {output_file}"
    assert Path(output_file).stat().st_size > 0, "Output file is empty"

    output_columns = get_mtz_columns(output_file)
    print(f"Output columns: {output_columns}")

    # Should have F+, SIGF+, F-, SIGF- (or Fplus, SIGFplus, Fminus, SIGFminus)
    expected = ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus']
    for col in expected:
        assert col in output_columns, f"Expected column {col} not found in {output_columns}"

    print("✅ Output contains FPAIR columns (Fplus, SIGFplus, Fminus, SIGFminus)")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
@pytest.mark.skipif(
    not check_ccp4_available(),
    reason="CCP4 not available. Run: source /Applications/ccp4-*/bin/ccp4.setup-sh"
)
def test_ipair_to_imean_conversion(tmp_path):
    """
    Test IPAIR → IMEAN conversion.

    Converts anomalous intensities (I+/I-) to mean intensities (I, SIGI)
    by averaging I+ and I-.
    """
    from ccp4i2.core.CCP4XtalData import CObsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Load IPAIR file (anomalous intensities)
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz"
    )
    assert Path(input_file).exists(), f"Input file not found: {input_file}"

    # Create file object
    ipair_file = CObsDataFile()
    ipair_file.setFullPath(input_file)
    ipair_file.setContentFlag()

    print(f"\nInput file: {input_file}")
    print(f"Input contentFlag: {ipair_file.contentFlag}")
    print(f"Input columns: {get_mtz_columns(input_file)}")

    assert int(ipair_file.contentFlag) == CObsDataFile.CONTENT_FLAG_IPAIR, \
        f"Expected IPAIR (1), got {ipair_file.contentFlag}"

    # Convert to IMEAN
    print("\n" + "=" * 60)
    print("Converting IPAIR → IMEAN...")
    print("=" * 60)

    output_file = ipair_file.as_IMEAN(work_directory=str(tmp_path))

    print(f"\n✅ Conversion completed: {output_file}")

    # Verify output exists and has correct columns
    assert Path(output_file).exists(), f"Output file not created: {output_file}"
    assert Path(output_file).stat().st_size > 0, "Output file is empty"

    output_columns = get_mtz_columns(output_file)
    print(f"Output columns: {output_columns}")

    # Should have I, SIGI
    assert 'I' in output_columns, f"Expected column I not found in {output_columns}"
    assert 'SIGI' in output_columns, f"Expected column SIGI not found in {output_columns}"

    print("✅ Output contains IMEAN columns (I, SIGI)")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
@pytest.mark.skipif(
    not check_ccp4_available(),
    reason="CCP4 not available. Run: source /Applications/ccp4-*/bin/ccp4.setup-sh"
)
def test_fpair_to_fmean_conversion(tmp_path):
    """
    Test FPAIR → FMEAN conversion.

    First converts IPAIR → FPAIR to get an FPAIR file,
    then converts FPAIR → FMEAN.
    """
    from ccp4i2.core.CCP4XtalData import CObsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Step 1: Create FPAIR file from IPAIR
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz"
    )

    ipair_file = CObsDataFile()
    ipair_file.setFullPath(input_file)
    ipair_file.setContentFlag()

    print(f"\nStep 1: Converting IPAIR → FPAIR...")
    fpair_path = ipair_file.as_FPAIR(work_directory=str(tmp_path))
    print(f"FPAIR file: {fpair_path}")

    # Step 2: Convert FPAIR → FMEAN
    fpair_file = CObsDataFile()
    fpair_file.setFullPath(fpair_path)
    fpair_file.setContentFlag()

    print(f"FPAIR contentFlag: {fpair_file.contentFlag}")
    assert int(fpair_file.contentFlag) == CObsDataFile.CONTENT_FLAG_FPAIR, \
        f"Expected FPAIR (2), got {fpair_file.contentFlag}"

    print("\n" + "=" * 60)
    print("Step 2: Converting FPAIR → FMEAN...")
    print("=" * 60)

    output_file = fpair_file.as_FMEAN(work_directory=str(tmp_path))

    print(f"\n✅ Conversion completed: {output_file}")

    # Verify output exists and has correct columns
    assert Path(output_file).exists(), f"Output file not created: {output_file}"
    assert Path(output_file).stat().st_size > 0, "Output file is empty"

    output_columns = get_mtz_columns(output_file)
    print(f"Output columns: {output_columns}")

    # Should have F, SIGF
    assert 'F' in output_columns, f"Expected column F not found in {output_columns}"
    assert 'SIGF' in output_columns, f"Expected column SIGF not found in {output_columns}"

    print("✅ Output contains FMEAN columns (F, SIGF)")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
@pytest.mark.skipif(
    not check_ccp4_available(),
    reason="CCP4 not available. Run: source /Applications/ccp4-*/bin/ccp4.setup-sh"
)
def test_conversion_chain_ipair_to_fmean(tmp_path):
    """
    Test full conversion chain: IPAIR → IMEAN → FMEAN.

    This tests that we can perform multiple conversions in sequence.
    """
    from ccp4i2.core.CCP4XtalData import CObsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Start with IPAIR
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "merged_intensities_native.mtz"
    )

    print("\n" + "=" * 60)
    print("Testing conversion chain: IPAIR → IMEAN → FMEAN")
    print("=" * 60)

    # Step 1: IPAIR → IMEAN
    print("\nStep 1: IPAIR → IMEAN...")
    ipair_file = CObsDataFile()
    ipair_file.setFullPath(input_file)
    ipair_file.setContentFlag()

    imean_path = ipair_file.as_IMEAN(work_directory=str(tmp_path))
    print(f"✅ IMEAN file: {imean_path}")
    print(f"   Columns: {get_mtz_columns(imean_path)}")

    # Step 2: IMEAN → FMEAN
    print("\nStep 2: IMEAN → FMEAN...")
    imean_file = CObsDataFile()
    imean_file.setFullPath(imean_path)
    imean_file.setContentFlag()

    assert int(imean_file.contentFlag) == CObsDataFile.CONTENT_FLAG_IMEAN, \
        f"Expected IMEAN (3), got {imean_file.contentFlag}"

    fmean_path = imean_file.as_FMEAN(work_directory=str(tmp_path))
    print(f"✅ FMEAN file: {fmean_path}")

    output_columns = get_mtz_columns(fmean_path)
    print(f"   Columns: {output_columns}")

    # Verify final output
    assert 'F' in output_columns, f"Expected column F not found"
    assert 'SIGF' in output_columns, f"Expected column SIGF not found"

    print("\n✅ Full conversion chain successful: IPAIR → IMEAN → FMEAN")


if __name__ == "__main__":
    # Allow running this test directly
    import sys
    sys.exit(pytest.main([__file__, "-v", "-s"]))
