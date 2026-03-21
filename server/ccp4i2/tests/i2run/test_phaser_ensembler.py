import gemmi
from .utils import demoData, i2run, download
from .urls import pdbe_pdb


def test_gamma_two_models(cif8xfm):
    """Test phaser.ensembler with two models."""
    model = demoData("gamma", "gamma_model.pdb")
    args = ["phaser_ensembler"]
    args += ["--XYZIN_LIST", model]
    args += ["--XYZIN_LIST", model]
    with i2run(args) as job:
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"
