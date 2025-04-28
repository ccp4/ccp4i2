"""
Copyright (C) 2010 University of York, Leiden University
"""

from ....script import crank2_script


class crank2_phdmmb(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_phdmmb'
  INTERRUPTABLE = True
  INTERRUPTLABEL = 'I am happy with the backbone trace.  Stop tracing after the current cycle!'
  SHORTTASKTITLE = ''
  TASKVERSION = 0.01

  out_params = ["XYZOUT","FPHOUT","FPHOUT_HL","XYZOUT_HAND2","FPHOUT_HAND2","FPHOUT_HL_HAND2","F_SIGFanom_OUT","F_SIGFanom_OUT2","F_SIGFanom_OUT3","F_SIGFanom_OUT4","F_SIGF_OUT","XYZOUT_SUBSTR"]
  perform = ["Hand1Score","Hand2Score","CC"]
