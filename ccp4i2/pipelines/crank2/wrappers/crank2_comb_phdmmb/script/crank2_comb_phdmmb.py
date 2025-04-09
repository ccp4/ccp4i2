"""
Copyright (C) 2010 University of York, Leiden University
"""

from ....script import crank2_script


class crank2_comb_phdmmb(crank2_script.crank2):

  TASKMODULE = 'test'
#  TASKTITLE = 'Crank2 DM'
  TASKNAME = 'crank2_comb_phdmmb'
  SHORTTASKTITLE = ''
  INTERRUPTABLE = True
  INTERRUPTLABEL = 'I am happy with the model.  Stop building after finishing the current building cycle!'
#  TASKCOMMAND = 'crank2.py'
  TASKVERSION = 0.01

  out_params = ["XYZOUT","XYZOUT_SUBSTR", "FPHOUT_2FOFC", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
  perform = ["RFactor","RFree"]
