from collections import Counter
import numpy as np
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
    counts = Counter(output.values())
    distinct = len(counts)
    assert distinct > 2, (
        f"only {distinct} distinct flag values -- the working set should be "
        "segmented so the free set can be re-partitioned later"
    )

    # and the segments are even. This is the check that caught the divergence
    # in the first place: it passed on the 21 even segments freerflag produces
    # and failed when the native implementation started returning a binary
    # column, whose working flag holds ~95% of the reflections.
    expected = 1.0 / distinct
    for flag, count in sorted(counts.items()):
        share = count / total
        assert 0.5 * expected <= share <= 2.0 * expected, (
            f"segment {flag} holds {share:.1%} of reflections; "
            f"{distinct} segments should hold about {expected:.1%} each"
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


def test_baz2b_sca():
    """Binary-free Scalepack import: read_scalepack -> gemmi split, no
    scalepack2mtz/cmtzsplit. The .sca carries cell/SG in its header and no
    FreeR, so import_merged generates a fresh free set."""
    sca = demoData(
        "baz2b",
        "BAZ2BA_x828.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8281_scaled.sca",
    )
    args = ["import_merged", "--HKLIN", sca]
    with i2run(args) as job:
        obs = gemmi.read_mtz_file(str(job / "OBSOUT.mtz"))
        labels = [c.label for c in obs.columns]
        # anomalous intensities -> I(+/-) pair (CObsDataFile content flag 1)
        assert "Iplus" in labels and "Iminus" in labels, f"OBSOUT columns {labels}"

        free_mtz = gemmi.read_mtz_file(str(job / "FREEOUT.mtz"))
        free_mtz.ensure_asu()
        freecol = free_mtz.rfree_column()
        assert freecol is not None, "FREEOUT.mtz missing FreeR column"
        flags = [int(f) for f in freecol]
        total = len(flags)
        free_fraction = sum(1 for f in flags if f == 0) / total
        assert 0.03 <= free_fraction <= 0.08, (
            f"generated test set is {free_fraction:.1%}; freerflag holds out ~5%"
        )
        assert len(set(flags)) > 2, "working set should be segmented"


# baz2b scalepack header: cell + C 2 2 21 (space group 20).
_BAZ2B_CELL = (83.09, 96.79, 57.95, 90.0, 90.0, 90.0)
_BAZ2B_SG = "C 2 2 21"


def _write_shelx_hklf4(sca_path, out_path):
    """Make a SHELX HKLF4 (intensities) file from a merged .sca, so we exercise
    the SHELX path without shipping a .hkl in the repo (all repo .hkl are XDS)."""
    from ccp4i2.lib.utils.files.reflection_formats import read_scalepack

    mtz = read_scalepack(sca_path, _BAZ2B_CELL, 20, anomalous=False)
    arr = np.array(mtz, copy=False)
    labels = [c.label for c in mtz.columns]
    hi, ki, li = labels.index("H"), labels.index("K"), labels.index("L")
    ii, si = labels.index("IMEAN"), labels.index("SIGIMEAN")
    with open(out_path, "w") as fh:
        for r in arr:
            I, S = r[ii], r[si]
            if I != I or S != S:   # skip missing (NaN)
                continue
            I = max(-9999.99, min(99999.99, float(I)))
            S = max(-9999.99, min(99999.99, float(S)))
            fh.write("%4d%4d%4d%8.2f%8.2f\n" % (int(r[hi]), int(r[ki]), int(r[li]), I, S))
        fh.write("%4d%4d%4d%8.2f%8.2f\n" % (0, 0, 0, 0, 0))


def test_shelx_hkl(tmp_path):
    """Binary-free SHELX import: read_shelx -> gemmi split, no f2mtz. SHELX
    carries neither cell/SG nor whether the columns are intensities or
    amplitudes, so all are supplied (cell + SG here, HKLF4 = intensities)."""
    sca = demoData(
        "baz2b",
        "BAZ2BA_x828.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8281_scaled.sca",
    )
    hkl = tmp_path / "baz2b_intensities.hkl"
    _write_shelx_hklf4(sca, hkl)

    a, b, c, al, be, ga = _BAZ2B_CELL
    args = [
        "import_merged", "--HKLIN", str(hkl),
        "--SPACEGROUP", _BAZ2B_SG,
        "--UNITCELL", f"a={a}", f"b={b}", f"c={c}",
        f"alpha={al}", f"beta={be}", f"gamma={ga}",
        "--SHELX_IS_INTENSITY", "True",
    ]
    with i2run(args) as job:
        obs = gemmi.read_mtz_file(str(job / "OBSOUT.mtz"))
        labels = [c.label for c in obs.columns]
        # intensities -> I/SIGI (mean), content flag 3
        assert "I" in labels and "SIGI" in labels, f"OBSOUT columns {labels}"

        free_mtz = gemmi.read_mtz_file(str(job / "FREEOUT.mtz"))
        free_mtz.ensure_asu()
        assert free_mtz.rfree_column() is not None, "FREEOUT.mtz missing FreeR"


def _write_merged_xds(sca_path, out_path):
    """Make a MERGED XDS_ASCII file from a merged .sca. The demo .hkl are all
    UNMERGED XDS (blocked), so a merged XDS to exercise the import path has to be
    generated. XDS carries its own cell/SG/wavelength."""
    from ccp4i2.lib.utils.files.reflection_formats import read_scalepack

    mtz = read_scalepack(sca_path, _BAZ2B_CELL, 20, anomalous=False)
    arr = np.array(mtz, copy=False)
    labels = [c.label for c in mtz.columns]
    hi, ki, li = labels.index("H"), labels.index("K"), labels.index("L")
    ii, si = labels.index("IMEAN"), labels.index("SIGIMEAN")
    a, b, c, al, be, ga = _BAZ2B_CELL
    with open(out_path, "w") as fh:
        fh.write("!FORMAT=XDS_ASCII    MERGE=TRUE    FRIEDEL'S_LAW=TRUE\n")
        fh.write("!SPACE_GROUP_NUMBER=20\n")
        fh.write(f"!UNIT_CELL_CONSTANTS= {a} {b} {c} {al} {be} {ga}\n")
        fh.write("!X-RAY_WAVELENGTH= 0.97950\n")
        fh.write("!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=5\n")
        fh.write("!ITEM_H=1\n!ITEM_K=2\n!ITEM_L=3\n!ITEM_IOBS=4\n!ITEM_SIGMA(IOBS)=5\n")
        fh.write("!END_OF_HEADER\n")
        for r in arr:
            I, S = r[ii], r[si]
            if I != I or S != S:
                continue
            fh.write("%6d%6d%6d %11.3E %11.3E\n" % (int(r[hi]), int(r[ki]), int(r[li]), I, S))
        fh.write("!END_OF_DATA\n")


def test_merged_xds():
    """Binary-free XDS import: gemmi read_xds_ascii -> to_mtz -> gemmi split.
    Merged XDS carries its own cell/SG, so nothing extra is supplied. (Unmerged
    XDS -- the usual CORRECT/INTEGRATE output -- is rejected by validity.)"""
    import tempfile
    import os

    sca = demoData(
        "baz2b",
        "BAZ2BA_x828.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8281_scaled.sca",
    )
    with tempfile.TemporaryDirectory() as d:
        xds = os.path.join(d, "merged_XDS_ASCII.HKL")
        _write_merged_xds(sca, xds)
        args = ["import_merged", "--HKLIN", xds]
        with i2run(args) as job:
            obs = gemmi.read_mtz_file(str(job / "OBSOUT.mtz"))
            labels = [c.label for c in obs.columns]
            assert "I" in labels and "SIGI" in labels, f"OBSOUT columns {labels}"
            free_mtz = gemmi.read_mtz_file(str(job / "FREEOUT.mtz"))
            free_mtz.ensure_asu()
            assert free_mtz.rfree_column() is not None, "FREEOUT.mtz missing FreeR"
