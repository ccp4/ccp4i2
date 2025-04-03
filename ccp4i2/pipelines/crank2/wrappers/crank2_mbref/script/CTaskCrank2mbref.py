"""
Copyright (C) 2011 University of York, Leiden University
"""

from ....script import CTaskCrank2


class CTaskCrank2mbref(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_mbref'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Model building'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 model building'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
