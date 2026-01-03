from ccp4i2.pipelines.crank2.script import crank2_script


class crank2_substrdet(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_substrdet'
  SHORTTASKTITLE = ''
  INTERRUPTABLE = True
  INTERRUPTLABEL = 'I am happy with the substructure.  Stop substructure detection!'
  TASKVERSION = 0.01

  out_params = ["XYZOUT_SUBSTR","XYZOUT_SUB_RES"]
  perform = ["CFOM",]
