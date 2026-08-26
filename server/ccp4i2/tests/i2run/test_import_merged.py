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
    """Check import_merged's outputs, and that the free set obeys freerflag's contract.

    Completing a free set is not a formatting question. What CCP4 freerflag
    guarantees, and what this asserts, is:

    - every reflection that was free stays free;
    - no reflection already used in refinement becomes free, which would
      quietly invalidate an R-free;
    - reflections that carried no flag are partitioned at the usual fraction;
    - the working set is segmented, so the set can be re-partitioned later.

    The same contract is checked directly against the CCP4 binary in
    tests/parity/test_freerflag_parity.py. This test previously asserted only
    that every flag value held between 2% and 15% of the reflections, which no
    binary free set can satisfy -- its working flag holds ~95% -- and which
    said nothing about the test set being preserved.
    """
    gemmi.read_mtz_file(str(job / "OBSOUT.mtz"))

    free_mtz = gemmi.read_mtz_file(str(job / "FREEOUT.mtz"))
    free_mtz.ensure_asu()
    freecol = free_mtz.rfree_column()
    assert freecol is not None, "FREEOUT.mtz missing FreeR column"

    hcol = free_mtz.column_with_label("H")
    kcol = free_mtz.column_with_label("K")
    lcol = free_mtz.column_with_label("L")
    output = {(h, k, l): int(f) for h, k, l, f in zip(hcol, kcol, lcol, freecol)}
    total = len(output)

    # the test set, at the fraction freerflag uses
    free_fraction = sum(1 for f in output.values() if f == 0) / total
    assert 0.03 <= free_fraction <= 0.08, (
        f"test set is {free_fraction:.1%} of reflections; freerflag holds out ~5%"
    )

    # the working set is segmented rather than left as a single flag
    distinct = len(set(output.values()))
    assert distinct > 2, (
        f"only {distinct} distinct flag values -- the working set should be "
        "segmented so the free set can be re-partitioned later"
    )

    # nothing crosses between the test and working sets
    input_flags = freer_flag_dict(freerin)
    was_free = [hkl for hkl, f in input_flags.items() if f == 0]
    was_working = [hkl for hkl, f in input_flags.items() if f != 0]
    if was_free:
        lost = [hkl for hkl in was_free if output.get(hkl, 0) != 0]
        assert not lost, (
            f"{len(lost)} of {len(was_free)} free reflections stopped being free"
        )
    if was_working:
        moved = [hkl for hkl in was_working if output.get(hkl, 1) == 0]
        assert not moved, (
            f"{len(moved)} reflections already used in refinement became free"
        )

    # COMPLETE mode extends the set, never shrinks it
    assert total >= len(input_flags), (
        f"output has fewer reflections ({total}) than input ({len(input_flags)})"
    )


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
