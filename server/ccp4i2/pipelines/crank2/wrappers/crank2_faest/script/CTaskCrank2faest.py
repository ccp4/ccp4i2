from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskCrank2faest(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_faest'
  TASKMODULE='test'
  TASKTITLE='FA estimation'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 FA estimation and other phasing preparations'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
