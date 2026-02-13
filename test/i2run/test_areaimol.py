from .urls import rcsb_pdb
from .utils import download, i2run
import gemmi

def test_from_pdb_rcsb():
    with download(rcsb_pdb("5a3h")) as pdb_5a3h:
        args = ["areaimol"]
        args += ["--XYZIN", pdb_5a3h]
        args += ["--OUTPUT_MODE", "ATOM"]
        with i2run(args) as job:
            mtz = gemmi.read_structure(str(job / "XYZOUT.pdb"))

def test_imol_from_pdb_rcsb():
    with download(rcsb_pdb("5a3h")) as pdb_5a3h:
        args = ["areaimol"]
        args += ["--XYZIN", pdb_5a3h]
        args += ["--DIFFMODE", "IMOL"]
        args += ["--OUTPUT_MODE", "ATOM"]
        with i2run(args) as job:
            mtz = gemmi.read_structure(str(job / "XYZOUT.pdb"))

def test_compare_from_pdb_rcsb():
    with download(rcsb_pdb("5a3h")) as pdb_5a3h, download(rcsb_pdb("8a3h")) as pdb_8a3h:
        args = ["areaimol"]
        args += ["--XYZIN", pdb_5a3h]
        args += ["--XYZIN2", pdb_8a3h]
        args += ["--DIFFMODE", "COMPARE"]
        args += ["--OUTPUT_MODE", "ATOM"]
        with i2run(args) as job:
            mtz = gemmi.read_structure(str(job / "XYZOUT.pdb"))
