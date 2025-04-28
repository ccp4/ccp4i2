"""
Copyright (C) 2013 STFC
"""

import glob
import os
import shutil

from ......core import CCP4PluginScript


class xia2_pointless(CCP4PluginScript.CPluginScript):

    TASKMODULE = None      # Where this plugin will appear on the gui
    TASKTITLE = 'Pointless in XIA2'
    TASKNAME = 'xia2_pointless'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin

    ERROR_CODES = { 101 : {'description' : 'XIA2 run directory does not exist' },
                    102 : {'description' : 'The XIA2 job failed with error file' },
                    103 : {'description' : 'No XIA2 run report file found' },
                    104 : {'description' : 'No XIA2 run output experimental data file found' }
                    }

    def process(self):
      xdir = self.container.inputData.XIA2_DIRECTORY.__str__()
      if not os.path.exists(xdir):
        self.errReport.append(self.__class__,101,xdir)
        self.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
        
      logfileList = glob.glob(os.path.join(xdir,'LogFiles','*pointless.log'))
      if len(logfileList)>0:
        shutil.copyfile(logfileList[0],os.path.join(self.workDirectory,'log.txt'))

      self.reportStatus(CCP4PluginScript.CPluginScript.SUCCEEDED)
