from shutil import which
import gemmi
from pytest import mark
from .utils import i2run


@mark.skipif(which("zanuda") is None and which("zanuda.exe") is None, reason="zanuda not installed")
def test_zanuda_8xfm(cif8xfm, mtz8xfm):
    """Test zanuda space group validation with 8xfm mmCIF input."""
    args = ["zanuda"]
    args += ["--XYZIN", f"fullPath={cif8xfm}"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]

    with i2run(args) as job:
        # Zanuda writes output as zanuda.pdb / zanuda.mtz (not XYZOUT)
        xyzout = job / "zanuda.pdb"
        assert xyzout.exists(), f"No zanuda.pdb: {list(job.iterdir())}"
        gemmi.read_pdb(str(xyzout))

        # Check split map coefficients
        for name in ["FPHIOUT", "DIFFPHIOUT"]:
            mtz_file = job / f"{name}.mtz"
            assert mtz_file.exists(), f"No {name}: {list(job.iterdir())}"

        # Each output says what it is and in which space group (they were
        # listed as 'zanuda.pdb', 'FPHIOUT.mtz' and 'DIFFPHIOUT.mtz').
        from ccp4i2.db import models
        job_record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
        for param in ("XYZOUT", "FPHIOUT", "DIFFPHIOUT"):
            gleaned = models.File.objects.filter(job=job_record, job_param_name=param).first()
            assert gleaned is not None and "Zanuda" in (gleaned.annotation or ""), \
                f"{param} has no annotation"
