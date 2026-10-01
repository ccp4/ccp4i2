import pytest
import gemmi
from .utils import demoData, i2run


def annotation(job, param):
    from ccp4i2.db import models
    record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
    gleaned = models.File.objects.filter(job=record, job_param_name=param).first()
    return gleaned.annotation if gleaned else None


def test_gamma_merge_two_mtz():
    """Test merging two mini-MTZ files from gamma data."""
    mtz1 = demoData("gamma", "merged_intensities_native.mtz")
    mtz2 = demoData("gamma", "merged_intensities_Xe.mtz")
    args = ["mergeMtz"]
    args += ["--MINIMTZINLIST", f"fileName={mtz1}"]
    args += ["--MINIMTZINLIST", f"fileName={mtz2}"]
    with i2run(args) as job:
        hklout = job / "HKLOUT.mtz"
        assert hklout.exists(), f"No HKLOUT: {list(job.iterdir())}"
        labels = [c.label for c in gemmi.read_mtz_file(str(hklout)).columns]
        # Both files are anomalous intensities. It wrote H,K,L alone (the
        # columns were asked for with the Qt-era columnNames signature); the
        # second file's names, already used, are prefixed by its position.
        assert labels == ["H", "K", "L", "Iplus", "SIGIplus", "Iminus", "SIGIminus",
                          "2_Iplus", "2_SIGIplus", "2_Iminus", "2_SIGIminus"], labels
        # It was listed as "HKLOUT.mtz".
        assert annotation(job, "HKLOUT") == ("Merged from 2 files: columns Iplus,SIGIplus,"
            "Iminus,SIGIminus,2_Iplus,2_SIGIplus,2_Iminus,2_SIGIminus")


def test_gamma_merge_with_tags():
    """A tag prefixes its file's columns, as its tooltip says."""
    args = ["mergeMtz"]
    args += ["--MINIMTZINLIST", f"fileName={demoData('gamma', 'merged_intensities_native.mtz')}",
             "columnTag=NAT"]
    args += ["--MINIMTZINLIST", f"fileName={demoData('gamma', 'merged_intensities_Xe.mtz')}",
             "columnTag=XE"]
    with i2run(args) as job:
        labels = [c.label for c in gemmi.read_mtz_file(str(job / "HKLOUT.mtz")).columns]
        assert labels[3:] == ["NAT_Iplus", "NAT_SIGIplus", "NAT_Iminus", "NAT_SIGIminus",
                              "XE_Iplus", "XE_SIGIplus", "XE_Iminus", "XE_SIGIminus"], labels
