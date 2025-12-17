from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskCrank2refatompick(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_refatompick'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='Iterative atom picking'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 substructure improvement'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
