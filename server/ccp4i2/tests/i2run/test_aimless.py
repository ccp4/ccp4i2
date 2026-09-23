from pathlib import Path
import xml.etree.ElementTree as ET
from gemmi import read_mtz_file
from pytest import approx
from .utils import demoData, i2run


def check_result(job: Path, spacegroup, resolution, rmeas):
    for name in [
        "FREERFLAG",
        "HKLOUT_0-observed_data_asIMEAN",
        "HKLOUT_0-observed_data",
        "HKLOUT_unmerged",
    ]:
        read_mtz_file(str(job / f"{name}.mtz"))
    tree = ET.parse(job / "program.xml")
    assert tree.find(".//Result/Dataset/SpacegroupName").text == spacegroup
    assert float(tree.find(".//ResolutionHigh/Overall").text) == approx(resolution)
    assert float(tree.find(".//Rmeas/Overall").text) == approx(rmeas)

    # COMPLETE_MTZ: reconstructed observed data + FreeR, tracked as a top-level
    # output so the Export MTZ button can serve it (found by the generic
    # fallback) and it survives cleanup / project export (#247).
    complete = job / "COMPLETE_MTZ.mtz"
    assert complete.is_file(), f"COMPLETE_MTZ.mtz missing from {job}"
    mtz = read_mtz_file(str(complete))
    labels = {c.label for c in mtz.columns}
    assert "FREER" in labels, (
        f"COMPLETE_MTZ should carry the FreeR column; got {sorted(labels)}"
    )


def free_fraction(job: Path) -> float:
    """The share of reflections flagged free in the generated FreeR set."""
    import numpy as np
    mtz = read_mtz_file(str(job / "FREERFLAG.mtz"))
    flags = np.array(mtz.column_with_label("FREER"), dtype=float)
    flags = flags[~np.isnan(flags)]
    return float((flags == 0).sum() / flags.size)


def test_gamma():
    mtz = demoData("gamma", "gamma_native.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]
    with i2run(args) as job:
        check_result(job, "P 21 21 21", 1.81, 0.061)
        # FREER_FRACTION left alone: the def.xml default of 0.05 applies.
        assert free_fraction(job) == approx(0.05, abs=0.01)


def subjob_plugins(job: Path) -> list[str]:
    """The task each sub-job ran, in sub-job order."""
    names = []
    for params in sorted(job.glob("job_*/params.xml"),
                         key=lambda p: int(p.parent.name.split("_")[1])):
        tree = ET.parse(params)
        element = tree.find(".//pluginName")
        names.append(element.text if element is not None else "?")
    return names


def test_mdm2_autocutoff_runs_aimless_twice_and_the_rest_once():
    """The ordinary AUTOCUTOFF path, pinned so nobody "fixes" it away.

    Two Aimless runs is the whole point: the first estimates a resolution
    limit from the CC-half analysis, the second applies it. Only the runs
    that feed that estimate are repeated -- ctruncate and freerflag consume
    the final Aimless output and happen once.
    """
    mtz = demoData("mdm2", "mdm2_unmerged.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}",
            "--AUTOCUTOFF", "True"]
    with i2run(args) as job:
        plugins = subjob_plugins(job)
        assert plugins.count("aimless") == 2, (
            f"the estimate-then-apply pair is the design: {plugins}")
        assert plugins.count("ctruncate") == 1, f"{plugins}"
        assert plugins.count("freerflag") == 1, f"{plugins}"


def test_gamma_autocutoff_stops_when_no_cutoff_is_needed():
    """When the data already reach the edge there is no second Aimless run
    -- "anygood" in process_post_aimless, and the pipeline finishes during
    the first one.

    It used to finish and then carry on. The pipeline reports its verdict
    from inside a chain of ordinary calls, all of which return, so finishing
    unwound back into the run loop, which ran Aimless, phaser_analysis,
    ctruncate and freerflag a second time and replaced the first run's
    outputs -- applying the cutoff that this path exists to say was not
    needed. Two of the three unmerged demo datasets take this path.

    SubstituteLigand sets AUTOCUTOFF on every aimless_pipe it runs.
    """
    mtz = demoData("gamma", "gamma_native.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}",
            "--AUTOCUTOFF", "True"]
    with i2run(args) as job:
        plugins = subjob_plugins(job)
        assert plugins.count("aimless") == 1, (
            f"no cutoff was needed, so there is nothing to apply: {plugins}")
        assert plugins.count("ctruncate") == 1, f"{plugins}"
        assert plugins.count("freerflag") == 1, f"{plugins}"
        # One dataset in, one observed-data file out.
        assert len(list(job.glob("HKLOUT_*-observed_data.mtz"))) == 1


def test_gamma_freer_fraction():
    # FREER_FRACTION is what the GUI edits; it must reach freerflag's FRAC.
    mtz = demoData("gamma", "gamma_native.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}",
            "--FREER_FRACTION", "0.1"]
    with i2run(args) as job:
        assert free_fraction(job) == approx(0.1, abs=0.01)


def test_mdm2():
    mtz = demoData("mdm2", "mdm2_unmerged.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]
    with i2run(args) as job:
        check_result(job, "P 61 2 2", 1.25, 0.068)
