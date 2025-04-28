"""
Copyright (C) 2011 University of York, Leiden University
"""

from ....script import CTaskCrank2


class CTaskCrank2phdmmb(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_phdmmb'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Density modification, poly-Ala tracing'
  SHORTTASKTITLE = None
  DESCRIPTION='SHELXE density modification and poly-Ala tracing'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
