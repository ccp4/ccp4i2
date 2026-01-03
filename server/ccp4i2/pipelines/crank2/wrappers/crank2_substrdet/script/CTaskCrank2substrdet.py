from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskCrank2substrdet(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_substrdet'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Substructure detection'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 substructure detection'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
