from ccp4i2.pipelines.crank2.script import crank2_script


class crank2_substrdet(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_substrdet'
  TASKVERSION = 0.01

  out_params = ["XYZOUT_SUBSTR","XYZOUT_SUB_RES"]
  perform = ["CFOM",]
