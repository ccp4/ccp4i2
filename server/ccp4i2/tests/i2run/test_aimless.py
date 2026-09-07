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
