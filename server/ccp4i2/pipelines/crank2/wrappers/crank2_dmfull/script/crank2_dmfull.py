from ccp4i2.pipelines.crank2.script import crank2_script


class crank2_dmfull(crank2_script.crank2):

  TASKMODULE = 'test'
#  TASKTITLE = 'Crank2 DM'
  TASKNAME = 'crank2_dmfull'
#  TASKCOMMAND = 'crank2.py'
  SHORTTASKTITLE = ''
  TASKVERSION = 0.01

  #out_params = ["XYZOUT_SUBSTR", "FPHOUT", "FPHOUT_DIFFANOM"]
  out_params = ["XYZOUT_SUBSTR", "FPHOUT", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
  perform = ["FOM",]
