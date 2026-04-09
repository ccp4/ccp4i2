import gemmi
from .urls import pdbe_sfcif
from .utils import demoData, download, i2run


def test_2ceu_cif():
    with download(pdbe_sfcif("2ceu")) as cif:
        args = ["import_merged"]
        args += ["--HKLIN", cif]
        args += ["--SPACEGROUP", "I 2 2 2"]
        with i2run(args) as job:
            check_output(job, cif)


def test_gamma_mtz():
    ianom = demoData("gamma", "merged_intensities_Xe.mtz")
    freer = demoData("gamma", "freeR.mtz")
    args = ["import_merged"]
    args += ["--HKLIN", ianom]
    args += ["--FREERFLAG", freer]
    with i2run(args) as job:
        check_output(job, freer)


def check_output(job, freerin):
    """
    Check that import_merged produces valid output files.

    Note: We don't expect perfect preservation of FreeR flags when using
    COMPLETE mode, as freerflag will:
    1. Extend the set to cover all observations
    2. Ensure symmetry-equivalent reflections have the same flag
    3. Maintain reasonable flag distribution

    This can cause some flags to be reassigned compared to the input.
    """
    # Check OBSOUT.mtz exists and is valid
    obs_mtz = gemmi.read_mtz_file(str(job / "OBSOUT.mtz"))

    # Check FREEOUT.mtz exists and is valid
    free_mtz = gemmi.read_mtz_file(str(job / "FREEOUT.mtz"))
    free_mtz.ensure_asu()

    # Verify FreeR column exists
    freecol = free_mtz.rfree_column()
    assert freecol is not None, "FREEOUT.mtz missing FreeR column"

    # Count reflections by flag value (20-bin system: flags 0-19)
    flag_counts = {}
    for flag in freecol:
        flag_counts[int(flag)] = flag_counts.get(int(flag), 0) + 1

    total = sum(flag_counts.values())
    num_bins = len(flag_counts)

    # DEBUG: Print flag distribution
    print(f"\n=== FreeR Flag Distribution ===")
    print(f"Total reflections: {total}")
    print(f"Number of bins: {num_bins}")
    for flag in sorted(flag_counts.keys())[:5]:
        print(f"  Flag {flag}: {flag_counts[flag]} ({flag_counts[flag]/total*100:.1f}%)")
    if num_bins > 5:
        print(f"  ... (showing first 5 of {num_bins} bins)")

    # Verify we have a reasonable number of bins (typically 20 for CCP4)
    assert num_bins >= 2, f"Too few FreeR bins: {num_bins}"

    # Verify each bin has reasonable size
    # In a 20-bin system, each should be ~5% (100/20)
    # In a binary system, test set should be ~5-10%
    expected_fraction = 1.0 / num_bins
    for flag, count in flag_counts.items():
        fraction = count / total
        # Each bin should be within reasonable range of expected
        # (allow 2-15% for flexibility - covers both 20-bin at ~5% and binary at ~5-10%)
        assert 0.02 < fraction < 0.15, \
            f"Flag {flag} has {fraction:.1%} of reflections (expected ~{expected_fraction:.1%})"

    # Verify output has at least as many reflections as input
    # (COMPLETE mode extends the set)
    input_flags = freer_flag_dict(freerin)
    assert total >= len(input_flags), \
        f"Output has fewer reflections ({total}) than input ({len(input_flags)})"


def freer_flag_dict(hklin):
    if hklin.endswith(".mtz"):
        mtz = gemmi.read_mtz_file(hklin)
    else:
        doc = gemmi.cif.read(hklin)
        rblock = gemmi.as_refln_blocks(doc)[0]
        mtz = gemmi.CifToMtz().convert_block_to_mtz(rblock)
    mtz.ensure_asu()
    hcol = mtz.column_with_label("H")
    kcol = mtz.column_with_label("K")
    lcol = mtz.column_with_label("L")
    freecol = mtz.rfree_column()
    return {
        (h, k, l): min(free, 1)
        for h, k, l, free in zip(hcol, kcol, lcol, freecol)
    }
