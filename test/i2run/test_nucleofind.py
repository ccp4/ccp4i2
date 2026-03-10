from gemmi import read_ccp4_map

from .urls import redo_mtz
from .utils import download, i2run


def test_1hr2():
    with download(redo_mtz("1hr2")) as mtz:
        args = ["nucleofind"]
        args += ["--FPHIIN", f"fullPath={mtz}", "columnLabels=/*/*/[FWT,PHWT]"]
        with i2run(args) as job:
            read_ccp4_map(str(job / "PHOSPHATE.map"))
            read_ccp4_map(str(job / "SUGAR.map"))
            read_ccp4_map(str(job / "BASE.map"))
