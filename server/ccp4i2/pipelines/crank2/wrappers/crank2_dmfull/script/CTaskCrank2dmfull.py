from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskCrank2dmfull(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_dmfull'
  TASKMODULE='test'
  TASKTITLE='Density modification'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 density modification'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
