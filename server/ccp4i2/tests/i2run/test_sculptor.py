from .utils import demoData, i2run


def test_gamma_sculptor():
    """Test sculptor model trimming with gamma alignment."""
    args = ["sculptor"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--ALIGNMENTORSEQUENCEIN", "ALIGNMENT"]
    args += ["--ALIGNIN", demoData("gamma", "gamma.pir")]
    with i2run(args) as job:
        # sculptor outputs a list of PDB files; check at least one exists
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"
