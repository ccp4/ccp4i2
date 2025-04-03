"""
Copyright (C) 2013 STFC
"""

import glob
import os
import shutil

from ......core import CCP4PluginScript


class xia2_integration(CCP4PluginScript.CPluginScript):

    TASKMODULE = None      # Where this plugin will appear on the gui
    TASKTITLE = 'Integration in XIA2'
    TASKNAME = 'xia2_integration'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin

    ERROR_CODES = { 101 : {'description' : 'XIA2 run directory does not exist' }                  
                    }

    def process(self):
      xdir = self.container.inputData.XIA2_DIRECTORY.__str__()
      if not os.path.exists(xdir):
        self.appendErrReport(cls=self.__class__,code=101,details=xdir)
        self.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
        
      datafileList = glob.glob(os.path.join(xdir,'DataFiles','Integrate','*.mtz'))
      datafileList.extend(glob.glob(os.path.join(xdir,'DataFiles','Integrate','*.sca')))
      datafileList.extend(glob.glob(os.path.join(xdir,'DataFiles','Integrate','*.HKL')))
      
      print('xia2_integration.process datafileList',datafileList)
      if len(datafileList)>0:
        fileName = os.path.join(self.workDirectory,'Unmerged'+os.path.splitext(datafileList[0])[1])
        shutil.copyfile(datafileList[0],fileName)
        self.container.outputData.HKLOUT.setFullPath(fileName)
        self.container.outputData.HKLOUT.annotation.set('Unmerged data')

      logfileList = glob.glob(os.path.join(xdir,'LogFiles','*integrate.log'))
      logfileList.extend(glob.glob(os.path.join(xdir,'LogFiles','*INTEGRATE.log')))
      if len(logfileList)>0:
        shutil.copyfile(logfileList[0],os.path.join(self.workDirectory,'log.txt'))
      self.reportStatus(CCP4PluginScript.CPluginScript.SUCCEEDED)
