import xml.etree.ElementTree as ET
from pathlib import Path
import gemmi
from .utils import demoData, i2run


def test_dimple():
    from ccp4i2.db import models

    args = ["i2Dimple"]  # Plugin is named i2Dimple in registry
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_native.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["FPHIOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")

        # COMPLETE_MTZ: dimple's unsplit reflection file (final.mtz), tracked as
        # a top-level output so the Export MTZ button can serve it (via the
        # generic fallback over tracked output-MTZ Files). A standalone wrapper
        # keeps the file's native name (final.mtz) rather than renaming to
        # COMPLETE_MTZ.mtz - the invariant is that it lives in the job's own
        # top-level dir, which it does. Assert on the gleaned File the export
        # path actually resolves, not a fixed filename.
        job_number = job.name.replace("job_", "")
        job_record = models.Job.objects.filter(number=job_number).first()
        complete = models.File.objects.filter(
            job=job_record, job_param_name="COMPLETE_MTZ"
        ).first()
        assert complete is not None, "COMPLETE_MTZ not gleaned as a tracked output"
        complete_path = Path(complete.path)
        assert complete_path.is_file(), (
            f"tracked COMPLETE_MTZ {complete_path} missing on disk"
        )
        assert complete_path.parent == job, (
            f"COMPLETE_MTZ must live in the job's own dir; got {complete_path}"
        )
        gemmi.read_mtz_file(str(complete_path))
