import json
import gemmi
import pytest
from .utils import download, i2run


_ROBETTA_MODELS = "https://robetta.bakerlab.org/models_download.php"


@pytest.fixture(scope="module", name="alphafold_info")
def alphafold_info_fixture():
    url = "https://alphafold.ebi.ac.uk/api/prediction/Q8W3K0"
    with download(url) as path:
        with open(path, encoding="utf-8") as f:
            yield json.load(f)[0]


@pytest.fixture(scope="module", name="alphafold_cif")
def alphafold_cif_fixture(alphafold_info):
    with download(alphafold_info["cifUrl"]) as path:
        yield path


@pytest.fixture(scope="module", name="alphafold_pdb")
def alphafold_pdb_fixture(alphafold_info):
    with download(alphafold_info["pdbUrl"]) as path:
        yield path


@pytest.fixture(scope="module", name="alphafold_pae")
def alphafold_pae_fixture(alphafold_info):
    with download(alphafold_info["paeDocUrl"]) as path:
        yield path


@pytest.fixture(scope="module", name="robetta_pdb")
def robetta_pdb_fixture():
    with download(f"{_ROBETTA_MODELS}?id=14697") as path:
        yield path


def test_alphafold_pdb(alphafold_pdb):
    args = ["editbfac"]
    args += ["--XYZIN", alphafold_pdb]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "converted_model.pdb"))
        gemmi.read_pdb(str(job / "converted_model_chainA1.pdb"))


def test_alphafold_cif(alphafold_cif):
    args = ["editbfac"]
    args += ["--XYZIN", alphafold_cif]
    args += ["--CONFCUT", "False"]
    args += ["--COMPACTREG", "False"]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "converted_model.pdb"))
        assert not (job / "converted_model_chainA1.pdb").exists()


def test_alphafold_pae(alphafold_cif, alphafold_pae):
    args = ["editbfac"]
    args += ["--XYZIN", alphafold_cif]
    args += ["--PAEIN", alphafold_pae]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "converted_model.pdb"))
        for i in range(1, 5):
            gemmi.read_pdb(str(job / f"converted_model_chainA{i}.pdb"))


def test_robetta(robetta_pdb):
    args = ["editbfac"]
    args += ["--XYZIN", robetta_pdb]
    args += ["--BTREATMENT", "rmsd"]
    args += ["--CONFCUT", "False"]
    args += ["--COMPACTREG", "False"]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "converted_model.pdb"))
        assert not (job / "converted_model_chainA1.pdb").exists()
