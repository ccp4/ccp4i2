from .utils import demoData, i2run


def test_gamma():
    """Test findmyseq sequence identification with gamma data."""
    args = ["findmyseq"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--F_SIGF", demoData("gamma", "gamma_native.mtz")]
    args += ["--FPHI", demoData("gamma", "initial_phases.mtz")]
    with i2run(args) as job:
        seqout = job / "SEQOUT.fasta"
        assert seqout.exists(), f"No SEQOUT: {list(job.iterdir())}"
        content = seqout.read_text()
        assert len(content) > 0, "Empty sequence output"
