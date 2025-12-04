"""
Test CPhsDataFile contentFlag conversions.

This test suite verifies all supported conversions between
MTZ phase data formats:
- HL (1): Hendrickson-Lattman coefficients (HLA, HLB, HLC, HLD)
- PHIFOM (2): Phase + Figure of Merit (PHI, FOM)

Conversion matrix:
                TO
          HL  PHIFOM
FROM HL    ✓     ✓
     PHIFOM ✓     ✓

Implementation uses gemmi's native HL coefficient support and
numerical methods for phase probability distributions.

Requires:
- CCP4I2_ROOT environment variable
- gemmi library (for MTZ I/O)
- numpy (for numerical calculations)
"""

import pytest
import os
import numpy as np
from pathlib import Path


def get_mtz_columns(mtz_path):
    """Get column names and types from an MTZ file using gemmi."""
    import gemmi
    mtz = gemmi.read_mtz_file(str(mtz_path))
    return [(col.label, col.type) for col in mtz.columns]


def get_mtz_data(mtz_path, column_label):
    """Extract a column's data from an MTZ file."""
    import gemmi
    mtz = gemmi.read_mtz_file(str(mtz_path))
    col = mtz.column_with_label(column_label)
    return np.array(col)


def check_fom_range(mtz_path):
    """Verify FOM values are in valid [0, 1] range."""
    import gemmi
    mtz = gemmi.read_mtz_file(str(mtz_path))

    # Find FOM column (type W)
    fom_col = None
    for col in mtz.columns:
        if col.type == 'W':
            fom_col = col
            break

    if fom_col is None:
        raise ValueError(f"No FOM column (type W) found in {mtz_path}")

    fom = np.array(list(fom_col), dtype=np.float32)
    valid = (fom >= 0.0) & (fom <= 1.0)
    return np.all(valid), fom.min(), fom.max()


def check_phi_range(mtz_path):
    """Verify PHI values are in valid [-180, 360] range (supports both conventions)."""
    import gemmi
    mtz = gemmi.read_mtz_file(str(mtz_path))

    # Find PHI column (type P)
    phi_col = None
    for col in mtz.columns:
        if col.type == 'P':
            phi_col = col
            break

    if phi_col is None:
        raise ValueError(f"No PHI column (type P) found in {mtz_path}")

    phi = np.array(list(phi_col), dtype=np.float32)
    # Accept both [0, 360] and [-180, 180] ranges
    valid = (phi >= -180.0) & (phi <= 360.0)
    return np.all(valid), phi.min(), phi.max()


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
def test_hl_to_phifom_conversion(tmp_path):
    """
    Test HL → PHIFOM conversion.

    Converts Hendrickson-Lattman coefficients (HLA, HLB, HLC, HLD) to
    best phase estimate (PHI) and figure of merit (FOM) using numerical
    optimization of the phase probability distribution.
    """
    from core.CCP4XtalData import CPhsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Load HL file (Hendrickson-Lattman coefficients)
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "initial_phases.mtz"
    )
    assert Path(input_file).exists(), f"Input file not found: {input_file}"

    # Create file object
    hl_file = CPhsDataFile()
    hl_file.setFullPath(input_file)
    hl_file.setContentFlag()

    print(f"\nInput file: {input_file}")
    print(f"Input contentFlag: {hl_file.contentFlag}")
    print(f"Input columns:")
    for label, col_type in get_mtz_columns(input_file):
        print(f"  {label:15s} type={col_type}")

    assert int(hl_file.contentFlag) == CPhsDataFile.CONTENT_FLAG_HL, \
        f"Expected HL (1), got {hl_file.contentFlag}"

    # Convert to PHIFOM
    print("\n" + "=" * 60)
    print("Converting HL → PHIFOM...")
    print("=" * 60)

    output_file = hl_file.as_PHIFOM(work_directory=str(tmp_path))

    print(f"\n✅ Conversion completed: {output_file}")

    # Verify output exists
    assert Path(output_file).exists(), f"Output file not created: {output_file}"
    assert Path(output_file).stat().st_size > 0, "Output file is empty"

    print(f"Output columns:")
    output_columns = get_mtz_columns(output_file)
    for label, col_type in output_columns:
        print(f"  {label:15s} type={col_type}")

    # Should have PHI (type P) and FOM (type W)
    # Note: chltofom may create columns with names like "PHI,FOM.Phi_fom.phi"
    column_labels = [label for label, _ in output_columns]
    col_types = {label: col_type for label, col_type in output_columns}

    # Find PHI column (type P)
    phi_cols = [label for label, ctype in output_columns if ctype == 'P']
    assert len(phi_cols) > 0, f"Expected phase column (type P) not found in {output_columns}"

    # Find FOM column (type W)
    fom_cols = [label for label, ctype in output_columns if ctype == 'W']
    assert len(fom_cols) > 0, f"Expected FOM column (type W) not found in {output_columns}"

    print(f"✅ Output contains PHIFOM columns:")
    print(f"   Phase (P): {phi_cols[0]}")
    print(f"   FOM (W):   {fom_cols[0]}")

    # Validate FOM range [0, 1]
    fom_valid, fom_min, fom_max = check_fom_range(output_file)
    assert fom_valid, f"FOM values out of range [0, 1]: min={fom_min}, max={fom_max}"
    print(f"✅ FOM values in valid range [0, 1]: min={fom_min:.3f}, max={fom_max:.3f}")

    # Validate PHI range (allow both [0,360] and [-180,180])
    phi_valid, phi_min, phi_max = check_phi_range(output_file)
    assert phi_valid, f"PHI values out of range [-180, 360]: min={phi_min}, max={phi_max}"
    print(f"✅ PHI values in valid range: min={phi_min:.1f}°, max={phi_max:.1f}°")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
def test_phifom_to_hl_conversion(tmp_path):
    """
    Test PHIFOM → HL conversion.

    First converts HL → PHIFOM to get a PHIFOM file,
    then converts PHIFOM → HL using centrosymmetric approximation.
    """
    from core.CCP4XtalData import CPhsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Step 1: Create PHIFOM file from HL
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "initial_phases.mtz"
    )

    hl_file = CPhsDataFile()
    hl_file.setFullPath(input_file)
    hl_file.setContentFlag()

    print(f"\nStep 1: Converting HL → PHIFOM...")
    phifom_path = hl_file.as_PHIFOM(work_directory=str(tmp_path))
    print(f"PHIFOM file: {phifom_path}")

    # Step 2: Convert PHIFOM → HL
    phifom_file = CPhsDataFile()
    phifom_file.setFullPath(phifom_path)
    phifom_file.setContentFlag()

    print(f"PHIFOM contentFlag: {phifom_file.contentFlag}")
    assert int(phifom_file.contentFlag) == CPhsDataFile.CONTENT_FLAG_PHIFOM, \
        f"Expected PHIFOM (2), got {phifom_file.contentFlag}"

    print("\n" + "=" * 60)
    print("Step 2: Converting PHIFOM → HL...")
    print("=" * 60)

    output_file = phifom_file.as_HL(work_directory=str(tmp_path))

    print(f"\n✅ Conversion completed: {output_file}")

    # Verify output exists
    assert Path(output_file).exists(), f"Output file not created: {output_file}"
    assert Path(output_file).stat().st_size > 0, "Output file is empty"

    print(f"Output columns:")
    output_columns = get_mtz_columns(output_file)
    for label, col_type in output_columns:
        print(f"  {label:15s} type={col_type}")

    # Should have HLA, HLB, HLC, HLD (type A)
    column_labels = [label for label, _ in output_columns]
    for hl_col in ['HLA', 'HLB', 'HLC', 'HLD']:
        assert hl_col in column_labels, f"Expected column {hl_col} not found"

    # Check column types
    col_types = {label: col_type for label, col_type in output_columns}
    for hl_col in ['HLA', 'HLB', 'HLC', 'HLD']:
        assert col_types[hl_col] == 'A', \
            f"{hl_col} column should be type A, got {col_types[hl_col]}"

    print("✅ Output contains HL columns (HLA, HLB, HLC, HLD all type=A)")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
def test_hl_phifom_roundtrip(tmp_path):
    """
    Test round-trip conversion: HL → PHIFOM → HL.

    This tests that we can perform multiple conversions in sequence
    and that the PHIFOM→HL conversion produces reasonable HL coefficients.
    """
    from core.CCP4XtalData import CPhsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    # Start with HL
    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "initial_phases.mtz"
    )

    print("\n" + "=" * 60)
    print("Testing round-trip: HL → PHIFOM → HL")
    print("=" * 60)

    # Step 1: HL → PHIFOM
    print("\nStep 1: HL → PHIFOM...")
    hl_file = CPhsDataFile()
    hl_file.setFullPath(input_file)
    hl_file.setContentFlag()

    # Get original HL data
    hla_orig = get_mtz_data(input_file, 'HLA')
    hlb_orig = get_mtz_data(input_file, 'HLB')

    phifom_path = hl_file.as_PHIFOM(work_directory=str(tmp_path))
    print(f"✅ PHIFOM file: {phifom_path}")

    # Step 2: PHIFOM → HL
    print("\nStep 2: PHIFOM → HL...")
    phifom_file = CPhsDataFile()
    phifom_file.setFullPath(phifom_path)
    phifom_file.setContentFlag()

    assert int(phifom_file.contentFlag) == CPhsDataFile.CONTENT_FLAG_PHIFOM, \
        f"Expected PHIFOM (2), got {phifom_file.contentFlag}"

    hl_reconstructed_path = phifom_file.as_HL(work_directory=str(tmp_path))
    print(f"✅ Reconstructed HL file: {hl_reconstructed_path}")

    # Verify reconstructed file has correct columns
    output_columns = get_mtz_columns(hl_reconstructed_path)
    column_labels = [label for label, _ in output_columns]

    for hl_col in ['HLA', 'HLB', 'HLC', 'HLD']:
        assert hl_col in column_labels, f"Missing column {hl_col}"

    print("✅ Reconstructed HL file has all HL columns")

    # Compare original vs reconstructed HL coefficients
    # Note: We use centrosymmetric approximation (HLC=HLD=0),
    # so we only compare HLA and HLB
    hla_recon = get_mtz_data(hl_reconstructed_path, 'HLA')
    hlb_recon = get_mtz_data(hl_reconstructed_path, 'HLB')
    hlc_recon = get_mtz_data(hl_reconstructed_path, 'HLC')
    hld_recon = get_mtz_data(hl_reconstructed_path, 'HLD')

    # HLC and HLD should be zero (centrosymmetric approximation)
    assert np.allclose(hlc_recon, 0.0, atol=1e-6), "HLC should be zero"
    assert np.allclose(hld_recon, 0.0, atol=1e-6), "HLD should be zero"
    print("✅ HLC and HLD are zero (as expected from centrosymmetric approximation)")

    # Calculate correlation between original and reconstructed HLA, HLB
    corr_a = np.corrcoef(hla_orig, hla_recon)[0, 1]
    corr_b = np.corrcoef(hlb_orig, hlb_recon)[0, 1]

    print(f"\nCorrelation coefficients (original vs reconstructed):")
    print(f"  HLA: {corr_a:.3f}")
    print(f"  HLB: {corr_b:.3f}")

    # We expect reasonable correlation since PHIFOM→HL uses an approximation
    # but may not be perfect due to information loss
    assert corr_a > 0.5, f"HLA correlation too low: {corr_a:.3f}"
    assert corr_b > 0.5, f"HLB correlation too low: {corr_b:.3f}"

    print("\n✅ Round-trip conversion successful with reasonable correlation")


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
def test_hl_numerical_accuracy(tmp_path):
    """
    Test numerical accuracy of HL → PHIFOM conversion.

    Verifies that the calculated PHI and FOM values are numerically
    stable and reasonable.
    """
    from core.CCP4XtalData import CPhsDataFile

    ccp4_root = os.environ["CCP4I2_ROOT"]

    input_file = os.path.join(
        ccp4_root, "demo_data", "gamma", "initial_phases.mtz"
    )

    hl_file = CPhsDataFile()
    hl_file.setFullPath(input_file)
    hl_file.setContentFlag()

    print("\n" + "=" * 60)
    print("Testing numerical accuracy of HL → PHIFOM")
    print("=" * 60)

    output_file = hl_file.as_PHIFOM(work_directory=str(tmp_path))

    # Extract converted data
    phi = get_mtz_data(output_file, 'PHI')
    fom = get_mtz_data(output_file, 'FOM')

    # Extract original HL coefficients
    hla = get_mtz_data(input_file, 'HLA')
    hlb = get_mtz_data(input_file, 'HLB')

    # Check for NaN or Inf values
    assert not np.any(np.isnan(phi)), "PHI contains NaN values"
    assert not np.any(np.isnan(fom)), "FOM contains NaN values"
    assert not np.any(np.isinf(phi)), "PHI contains Inf values"
    assert not np.any(np.isinf(fom)), "FOM contains Inf values"
    print("✅ No NaN or Inf values in output")

    # Check that FOM is higher for stronger phase information
    # When HLA and HLB are both large, we expect higher FOM
    hl_magnitude = np.sqrt(hla**2 + hlb**2)

    # Find reflections with strong vs weak phase information
    strong_idx = hl_magnitude > np.percentile(hl_magnitude, 75)
    weak_idx = hl_magnitude < np.percentile(hl_magnitude, 25)

    fom_strong = fom[strong_idx].mean()
    fom_weak = fom[weak_idx].mean()

    print(f"\nFOM statistics by HL magnitude:")
    print(f"  Strong phase info (top 25%): FOM mean = {fom_strong:.3f}")
    print(f"  Weak phase info (bottom 25%): FOM mean = {fom_weak:.3f}")

    # FOM should be higher for stronger phase information
    assert fom_strong > fom_weak, \
        f"Expected higher FOM for strong phases, got {fom_strong:.3f} vs {fom_weak:.3f}"

    print("✅ FOM correlates with HL coefficient magnitude (strong phases → higher FOM)")

    # Check phase relationship with HLA/HLB
    # PHI should be approximately atan2(HLB, HLA) for centrosymmetric approximation
    phi_approx = np.degrees(np.arctan2(hlb, hla))
    phi_approx = (phi_approx + 360) % 360  # Ensure [0, 360]

    # Calculate circular correlation (accounting for 360° periodicity)
    phase_diff = np.abs(phi - phi_approx)
    phase_diff = np.minimum(phase_diff, 360 - phase_diff)  # Handle wrap-around
    mean_phase_error = phase_diff.mean()

    print(f"\nPhase comparison with atan2(HLB, HLA) approximation:")
    print(f"  Mean absolute difference: {mean_phase_error:.1f}°")

    # Mean error should be reasonable (phases aren't exact due to HLC/HLD contributions)
    assert mean_phase_error < 45.0, \
        f"Mean phase error too large: {mean_phase_error:.1f}°"

    print("✅ Calculated phases reasonably match atan2(HLB, HLA)")


if __name__ == "__main__":
    # Allow running this test directly
    import sys
    sys.exit(pytest.main([__file__, "-v", "-s"]))
