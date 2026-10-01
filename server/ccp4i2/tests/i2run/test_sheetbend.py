import xml.etree.ElementTree as ET
from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm, mtz8xfm):
    from ccp4i2.db import models

    args = ["sheetbend"]
    args += ["--XYZIN", cif8xfm]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
        xml = ET.parse(job / "program.xml")
        rwork = float(xml.find(".//Final/Rwork").text)
        rfree = float(xml.find(".//Final/Rfree").text)
        assert rwork < 0.23
        assert rfree < 0.25
        # The model says what it is (it was listed as 'XYZOUT.cif').
        job_record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
        xyzout = models.File.objects.filter(job=job_record, job_param_name="XYZOUT").first()
        assert xyzout is not None and xyzout.annotation, "XYZOUT has no annotation"
