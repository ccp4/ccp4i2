from ccp4i2.pipelines.crank2.script import crank2_script


class crank2_comb_phdmmb(crank2_script.crank2):

  TASKNAME = 'crank2_comb_phdmmb'

  out_params = ["XYZOUT","XYZOUT_SUBSTR", "FPHOUT_2FOFC", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
  perform = ["RFactor","RFree"]
