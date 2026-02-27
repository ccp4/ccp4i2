from .utils import i2run


def test_8xfm(mtz8xfm, seq8xfm):
    args = ["AMPLE"]
    args += ["--AMPLE_F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--AMPLE_SEQIN", seq8xfm]
    args += ["--AMPLE_EXISTING_MODELS", "False"]
    args += ["--AMPLE_NPROC", "8"]
    i2run(args)
