"""
Copyright (C) 2011 University of York, Leiden University
"""

from ....script import CTaskCrank2


class CTaskCrank2dmfull(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_dmfull'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Density modification'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 density modification'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
