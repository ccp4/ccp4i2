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


def test_gamma():
    mtz = demoData("gamma", "gamma_native.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]
    with i2run(args) as job:
        check_result(job, "P 21 21 21", 1.81, 0.061)


def test_mdm2():
    mtz = demoData("mdm2", "mdm2_unmerged.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]
    with i2run(args) as job:
        check_result(job, "P 61 2 2", 1.25, 0.068)
