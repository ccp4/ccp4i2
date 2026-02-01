"""
     import_xia2.py: CCP4 GUI Project
     Copyright (C) 2013 STFC
"""

import os
from ccp4i2.core import CCP4PluginScript
from ccp4i2.pipelines.import_xia2.wrappers.xia2_run.script import xia2_run

class import_xia2(CCP4PluginScript.CPluginScript):

    TASKTITLE = 'Import XIA2 results'
    TASKNAME = 'import_xia2'
    WHATNEXT = ['aimless_pipe','phaser_pipeline','molrep_mr']

    ERROR_CODES = { 201 : {'description' : 'XIA2 directory does not exist' },
                    202 : {'description' : 'The XIA2 job failed with error file' },
                    203 : {'description' : 'Failed to find XIA2 run directories - is the chosen XIA2 directory correct?' }
                    }

    def process(self):
      xdir = self.container.inputData.XIA2_DIRECTORY.__str__()
      runModeList = []
      #print 'import_xia2.process',self.container.inputData.XIA2_DIRECTORY,xdir
      if not os.path.exists(xdir):
        self.appendErrorReport(201,xdir)
        self.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
        return CCP4PluginScript.CPluginScript.FAILED
      if os.path.exists(os.path.join(xdir,'xia2.error')):
        self.appendErrorReport(202,os.path.join(xdir,'xia2.error'))
        self.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
        return CCP4PluginScript.CPluginScript.FAILED
    
      if os.path.split(xdir)[1].endswith('-run'):
        runModeList = [ os.path.split(xdir)[1][0:-4] ]
        xdir = os.path.split(xdir)[0]
      else:
        import glob
        dirList = glob.glob(os.path.join(xdir,'*-run'))
        for d in dirList: runModeList.append(os.path.split(d)[1][0:-4] )
        if len(runModeList)==0:
          dirList = glob.glob(os.path.join(xdir,'2d*'))
          dirList.extend( glob.glob(os.path.join(xdir,'3d*')) )
          for d in dirList: runModeList.append(os.path.split(d)[1])
          

      #print 'runModeList',runModeList
      if len(runModeList)==0:
        self.appendErrorReport(203,xdir)
        return CCP4PluginScript.CPluginScript.FAILED
     
      makePluginObjectMode = 2
      for runMode in runModeList:
          if makePluginObjectMode == 2:
            pluginObj = self.makePluginObject('xia2_run',mode=makePluginObjectMode,pluginTitle='Import XIA2')
          else:
            pluginObj = self.makePluginObject('xia2_run',mode=makePluginObjectMode)
          #print 'pluginName',pluginName,pluginObj
          pluginObj.container.inputData.XIA2_DIRECTORY.setFullPath(xdir)
          pluginObj.container.controlParameters.RUN_MODE = runMode
          pluginObj.container.header.pluginTitle = pluginObj.runTitle(runMode)
          pluginObj.process()
          makePluginObjectMode = 1
          

      self.finished.emit(CCP4PluginScript.CPluginScript.SUCCEEDED)
      return CCP4PluginScript.CPluginScript.SUCCEEDED
