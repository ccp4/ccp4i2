import gemmi
from .utils import demoData, i2run


def test_gamma_model():
    """Test pdbset with a simple CRYST1 keyword."""
    args = ["pdbset_ui"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--EXTRA_PDBSET_KEYWORDS", "CHAIN A"]
    with i2run(args) as job:
        xyzout = job / "XYZOUT.pdb"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
        st = gemmi.read_structure(str(xyzout))
        assert len(st[0]) > 0
        # The gleaner records what the model is when the wrapper does not:
        # a PDB-format model was recorded with content 0 (not recognised).
        from ccp4i2.db import models
        record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
        out = models.File.objects.filter(job=record, job_param_name="XYZOUT").first()
        assert out.content == 1, out.content
        assert out.annotation == "Model edited by pdbset: CHAIN A", out.annotation
