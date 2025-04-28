"""
Copyright (C) 2011 University of York, Leiden University
"""

from ....script import CTaskCrank2


class CTaskCrank2faest(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_faest'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='FA estimation'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 FA estimation and other phasing preparations'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
