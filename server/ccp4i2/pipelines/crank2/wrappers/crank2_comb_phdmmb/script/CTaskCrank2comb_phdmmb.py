from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskCrank2comb_phdmmb(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_comb_phdmmb'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Model building'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 model building'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
