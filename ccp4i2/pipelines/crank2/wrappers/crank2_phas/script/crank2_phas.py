"""
Copyright (C) 2010 University of York, Leiden University
"""

from ....script import crank2_script


class crank2_phas(crank2_script.crank2):

  TASKMODULE = 'test'
#  TASKTITLE = 'Crank2 substrdetect'
  TASKNAME = 'crank2_phas'
  SHORTTASKTITLE = ''
  TASKVERSION = 0.01

  out_params = ["XYZOUT_SUBSTR", "FPHOUT", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
  perform = ["FOM",]
