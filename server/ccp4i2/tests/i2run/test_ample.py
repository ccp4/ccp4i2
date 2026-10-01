"""AMPLE with its helical ensembles (no models of one's own).

This test called i2run(args) as a statement, so the job never ran and it
"passed" in a tenth of a second. Run for real, AMPLE prepares its ensembles
and every MrBUMP search fails at once: ample/util/ample_util.py computes
ncopies = nresidues / nres, a float under Python 3 (15.7 here), and MrBUMP
rejects "LOCALFILE ... COPIES 15.7" (an integer is required). AMPLE then ends
with an empty summary and no model. An upstream bug in the CCP4 bundle
(ccp4-20260702 and ccp4-20260904); strict, so a fixed AMPLE shows up here.
"""
import pytest

from .utils import i2run


@pytest.mark.xfail(strict=True, reason="AMPLE writes a non-integer MrBUMP COPIES (Python 3 division)")
def test_8xfm(mtz8xfm, seq8xfm):
    args = ["AMPLE"]
    args += ["--AMPLE_F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--AMPLE_SEQIN", seq8xfm]
    args += ["--AMPLE_EXISTING_MODELS", "False"]
    args += ["--AMPLE_NPROC", "8"]
    with i2run(args) as job:
        assert list(job.glob("*.pdb")) or list(job.glob("*.cif")), list(job.iterdir())
