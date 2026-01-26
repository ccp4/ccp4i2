from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskCrank2phdmmb(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_phdmmb'
  TASKMODULE='test'
  TASKTITLE='Density modification, poly-Ala tracing'
  SHORTTASKTITLE = None
  DESCRIPTION='SHELXE density modification and poly-Ala tracing'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
