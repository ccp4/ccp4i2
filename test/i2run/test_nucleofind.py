from gemmi import read_ccp4_map
from pytest import fixture

from .urls import redo_mtz
from .utils import download, i2run


@fixture(scope="module", name="mtz")
def mtz_fixture():
    with download(redo_mtz("1hr2")) as path:
        yield path


def test_1hr2(mtz):
    args = ["nucleofind"]
    args += ["--FPHIIN", f"fullPath={mtz}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        read_ccp4_map(str(job / "PHOSPHATE.map"))
        read_ccp4_map(str(job / "SUGAR.map"))
        read_ccp4_map(str(job / "BASE.map"))


def test_1hr2_raw(mtz):
    args = ["nucleofind"]
    args += ["--FPHIIN", f"fullPath={mtz}", "columnLabels=/*/*/[FWT,PHWT]"]
    args += ["--OUTPUT_TYPE", "RAW"]
    with i2run(args) as job:
        read_ccp4_map(str(job / "PHOSPHATE.map"))
        read_ccp4_map(str(job / "SUGAR.map"))
        read_ccp4_map(str(job / "BASE.map"))
