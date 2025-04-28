"""
Copyright (C) 2010 University of York, Leiden University
"""

from ....script import crank2_script


class crank2_ref(crank2_script.crank2):

  TASKMODULE = 'test'
#  TASKTITLE = 'Crank2 DM'
  TASKNAME = 'crank2_ref'
  SHORTTASKTITLE = ''
#  TASKCOMMAND = 'crank2.py'
  TASKVERSION = 0.01

  out_params = ["XYZOUT","XYZOUT_SUBSTR", "FPHOUT_2FOFC", "FPHOUT_HL", "FPHOUT_DIFFANOM", "FPHOUT_DIFF"]
  perform = ["RFactor","RFree"]
