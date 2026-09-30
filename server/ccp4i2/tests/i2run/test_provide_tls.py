from .utils import demoData, i2run


def test_auto_from_coords():
    """Test ProvideTLS generating TLS groups from coordinates."""
    from ccp4i2.db import models

    args = ["ProvideTLS"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        tlsout = job / "TLSFILE.tls"
        assert tlsout.exists(), f"No TLS file: {list(job.iterdir())}"
        content = tlsout.read_text()
        assert "TLS" in content, "TLS file has no TLS keyword"
        # The file says what it holds (it was listed as 'TLSFILE.tls').
        job_record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
        tls = models.File.objects.filter(job=job_record, job_param_name="TLSFILE").first()
        assert tls is not None and tls.annotation.startswith("TLS definitions, 1 group")
