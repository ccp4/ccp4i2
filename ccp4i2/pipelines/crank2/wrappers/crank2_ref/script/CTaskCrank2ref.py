"""
Copyright (C) 2011 University of York, Leiden University
"""

from ....script import CTaskCrank2


class CTaskCrank2ref(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_ref'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Refinement'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 refinement'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
