"""
Copyright (C) 2010 University of York, Leiden University
"""

from ....script import crank2_script


class crank2_handdet(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_handdet'
  SHORTTASKTITLE = ''
  TASKVERSION = 0.01

  out_params = ["XYZOUT_SUBSTR","FPHOUT_HL","XYZOUT_HAND2","FPHOUT_HL_HAND2","F_SIGFanom_OUT","F_SIGFanom_OUT2","F_SIGFanom_OUT3","F_SIGFanom_OUT4","F_SIGF_OUT"]
  perform = ["Hand1Score","Hand2Score",]
