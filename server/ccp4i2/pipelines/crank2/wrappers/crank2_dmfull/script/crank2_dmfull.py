from ccp4i2.pipelines.crank2.script import crank2_script


class crank2_dmfull(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_dmfull'

  out_params = ["XYZOUT_SUBSTR", "FPHOUT", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
  perform = ["FOM",]
