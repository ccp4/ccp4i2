"""
Copyright (C) 2011 University of York, Leiden University
"""

from ....script import CTaskCrank2


class CTaskCrank2handdet(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_handdet'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Hand determination'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 determination of handedness'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
