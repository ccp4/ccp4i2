# Copyright (C) 2025 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
