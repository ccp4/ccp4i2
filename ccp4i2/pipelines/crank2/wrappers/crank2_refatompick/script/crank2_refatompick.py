"""
Copyright (C) 2010 University of York, Leiden University
"""

from ....script import crank2_script


class crank2_refatompick(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_refatompick'
  TASKVERSION = 0.01
  SHORTTASKTITLE = ''

  out_params = ["XYZOUT_SUBSTR", "XYZOUT", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
