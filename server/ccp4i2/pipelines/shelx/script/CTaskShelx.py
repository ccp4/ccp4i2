from ccp4i2.pipelines.crank2.script import CTaskCrank2


class CTaskShelx(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'shelx'
  TASKVERSION = 0.0
  TASKMODULE='expt_phasing'
  TASKTITLE='Automated structure solution - SHELXC/D/E phasing and building'
  SHORTTASKTITLE='SHELX'
  DESCRIPTION='Experimental phasing pipeline SHELX (run via Crank2)'
  RANK=1

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
