from ccp4i2.pipelines.crank2.script import crank2_script


class crank2_ref(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_ref'

  out_params = ["XYZOUT","XYZOUT_SUBSTR", "FPHOUT_2FOFC", "FPHOUT_HL", "FPHOUT_DIFFANOM", "FPHOUT_DIFF"]
  perform = ["RFactor","RFree"]
