from ccp4i2.pipelines.crank2.script import crank2_script


class crank2_mbref(crank2_script.crank2):

  TASKMODULE = 'test'
#  TASKTITLE = 'Crank2 DM'
  TASKNAME = 'crank2_mbref'
#  TASKCOMMAND = 'crank2.py'
  INTERRUPTABLE = True
  INTERRUPTLABEL = 'I am happy with the model.  Stop building after the current cycle!'
  TASKVERSION = 0.01
  SHORTTASKTITLE = ''

  out_params = ["XYZOUT","XYZOUT_SUBSTR", "FPHOUT_2FOFC", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
  perform = ["RFactor","RFree"]
