from pathlib import Path
from gemmi import read_structure
from pytest import fixture
import json
from .urls import pdbe_mmcif
from .utils import download, i2run


# TODO: Test long ligand names


@fixture(scope="module", name="cif")
def cif_fixture():
    with download(pdbe_mmcif("4dl8")) as path:
        yield path


def test_metalcoord(cif):
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
        with open(job / "AF3.json") as f:
            data = json.load(f)
        assert len(data) == 1
        assert len(data[0]["ligands"][0]["base"]) == 3
        assert len(data[0]["ligands"][0]["pdb"]) == 2      #  3 using default parameters 
        assert len(data[0]["ligands"][0]["angles"]) == 10  # 15 using default parameters
        assert data[0]["ligands"][0]["base"][0]["distance"][0] > 1.0
        assert data[0]["ligands"][0]["base"][0]["std"][0] > 0.01
        assert data[0]["ligands"][0]["base"][0].get("distance_model", None)
        assert data[0]["ligands"][0]["pdb"][0]["distance"][0] > 1.0
        assert data[0]["ligands"][0]["pdb"][0]["std"][0] > 0.01
        assert data[0]["ligands"][0]["pdb"][0].get("distance_model", None)
        assert data[0]["ligands"][0]["angles"][0].get("angle", None)
        assert data[0]["ligands"][0]["angles"][0].get("std", None)
        assert data[0]["ligands"][0]["angles"][0].get("angle_model", None)
        for ext in ("_coot.txt", ".mmcif", ".params", ".pdb", ".txt"):
            assert Path(job / f"AF3_restraints{ext}").exists()
        read_structure(str(job / "AF3_restraints.mmcif"))
        read_structure(str(job / "AF3_restraints.pdb"))
