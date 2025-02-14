"""
     crank2.py: CCP4 GUI Project
, Leiden University


"""

from ....script import crank2_script


class crank2_faest(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_faest'
  SHORTTASKTITLE = ''
  TASKVERSION = 0.01

  out_params = []
  perform = []

