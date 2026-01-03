from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskCrank2phas(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_phas'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Phasing'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 phasing'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
