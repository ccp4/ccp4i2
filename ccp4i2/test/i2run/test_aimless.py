from gemmi import read_mtz_file
from pytest import mark
from .utils import demoData, i2run


@mark.parametrize(
    "file",
    [
        demoData("gamma", "gamma_native.mtz"),
        demoData("mdm2", "mdm2_unmerged.mtz"),
    ],
)
def test_aimless(file: str):
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={file}"]
    with i2run(args) as job:
        for name in [
            "FREERFLAG",
            "HKLOUT_0-observed_data_asIMEAN",
            "HKLOUT_0-observed_data",
            "HKLOUT_unmerged",
        ]:
            read_mtz_file(str(job / f"{name}.mtz"))
