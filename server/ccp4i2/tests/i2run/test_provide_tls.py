from .utils import demoData, i2run


def test_auto_from_coords():
    """Test ProvideTLS generating TLS groups from coordinates."""
    args = ["ProvideTLS"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        tlsout = job / "TLSFILE.tls"
        assert tlsout.exists(), f"No TLS file: {list(job.iterdir())}"
        content = tlsout.read_text()
        assert "TLS" in content, "TLS file has no TLS keyword"
