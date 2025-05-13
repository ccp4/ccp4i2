from pathlib import Path
from gemmi import read_structure
from pytest import fixture
from .urls import pdbe_mmcif
from .utils import download, i2run


# TODO: Test long ligand names


@fixture(scope="module", name="cif")
def cif_fixture():
    with download(pdbe_mmcif("4dl8")) as path:
        yield path


def test_single_atom_mr(cif):
    args = ["metalCoord"]
    args += ["--XYZIN", cif]
    args += ["--LIGAND_CODE", "AF3"]
    args += ["--KEEP_LINKS", "True"]
    args += ["--MAXIMUM_COORDINATION_NUMBER", "5"]
    args += ["--COORD05", "trigonal-bipyramid"]
    args += ["--MINIMUM_SAMPLE_SIZE", "10"]
    args += ["--DISTANCE_THRESHOLD", "0.45"]
    args += ["--PROCRUSTES_DISTANCE_THRESHOLD", "0.2"]
    with i2run(args) as job:
        for ext in ("_coot.txt", ".mmcif", ".params", ".pdb", ".txt"):
            assert Path(job / f"AF3_restraints{ext}").exists()
        read_structure(str(job / "AF3_restraints.mmcif"))
        read_structure(str(job / "AF3_restraints.pdb"))
