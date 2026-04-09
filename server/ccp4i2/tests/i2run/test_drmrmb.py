# TODO: Test long ligand names (e.g. 8xfm)


# def test_drmrmb(mmcif, mtz, fasta):
#     args = ["dr_mr_modelbuild_pipeline"]
#     args += ["--XYZIN", mmcif]
#     args += ["--F_SIGF_IN", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREER_IN", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
#     args += ["--ASUIN", f"seqFile={fasta}"]
#     args += ["--MODELCRAFT_NCYC", "5"]
#     args += ["--REFMAC_NCYC", "0"]
#     args += ["--MERGED_OR_UNMERGED", "MERGED"]
#     args += ["--NMON", "1"]
