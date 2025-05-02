from pathlib import Path
import xml.etree.ElementTree as ET
from gemmi import read_mtz_file
from pytest import approx
from .utils import demoData, i2run


def check_result(job: Path):
    for name in [
        "FREERFLAG",
        "HKLOUT_0-observed_data_asIMEAN",
        "HKLOUT_0-observed_data",
        "HKLOUT_unmerged",
    ]:
        read_mtz_file(str(job / f"{name}.mtz"))
    tree = ET.parse(job / "XMLOUT.xml")
    return {
        "spacegroup": tree.find(".//Result/Dataset/SpacegroupName").text,
        "resolution": float(tree.find(".//ResolutionHigh/Overall").text),
        "rmeas": float(tree.find(".//Rmeas/Overall").text),
    }


def test_gamma():
    mtz = demoData("gamma", "gamma_native.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]
    with i2run(args) as job:
        result = check_result(job)
        assert result["spacegroup"] == "P 21 21 21"
        assert result["resolution"] == approx(1.81)
        assert result["rmeas"] == approx(0.061)


def test_mdm2():
    mtz = demoData("mdm2", "mdm2_unmerged.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]
    with i2run(args) as job:
        result = check_result(job)
        assert result["spacegroup"] == "P 61 2 2"
        assert result["resolution"] == approx(1.25)
        assert result["rmeas"] == approx(0.068)
