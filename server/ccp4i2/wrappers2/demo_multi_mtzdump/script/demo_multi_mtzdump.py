from __future__ import print_function

"""
    demo_multi_mtzdump.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

# Simple example of running external processes in blocking mode
# Multiple calls to mtzdump to find the cell in a list of mtzfiles

from ccp4i2.core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER # Utility to access CCP4i modules
from ccp4i2.core.CCP4PluginScript import CPluginScript

from ccp4i2.baselayer import QtCore
     
class demo_multi_mtzdump(CPluginScript):

    TASKMODULE = 'demo'
    TASKTITLE = 'Demo multi MTZ dump'
    TASKNAME = 'demo_multi_mtzdump'
    TASKVERSION= 0.0
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 0.1

    def process(self):
      # Run mtzdump on multiple mtzs to get the cell parameters
      import glob,os,time
      from ccp4i2.core.CCP4Utils import getCCP4I2Dir
      # Get a list on MTZ files to work on
      self.mtzlist = glob.glob(os.path.join(getCCP4I2Dir(),'wrappers','*','test_data','*.mtz'))
      #self.mtzlist.extend(glob.glob(os.path.join(getCCP4I2Dir(),'pipelines','*','test_data','*.mtz')))
      self.subProcessList = []
      self.startTime = time.time()

      # Loop over the MTZ files calling mtzdump for each one
      for mtz in self.mtzlist:
        # Create an instance of mtzdump plugin and set the HKLIN
        mtzdump = self.makePluginObject(pluginName='mtzdump',reportToDatabase=False)
        mtzdump.container.inputData.HKLIN.fullPath  = mtz
        # Set to run asynchronously and set a callback
        mtzdump.doAsync = True
        mtzdump.finished.connect(self.handleDone)
        #Start process
        rv = mtzdump.process()
        #The mtzdump instance must be saved to keep it in scope and everything else can be got from that.
        self.subProcessList.append(mtzdump)

      # Set a callback to end the entire process after a given time 
      self.timedCallback(demo_multi_mtzdump.TIMEOUT_PERIOD,self.handleTimeout)

      return CPluginScript.SUCCEEDED
        

    @QtCore.Slot()
    def handleDone(self,ret):
      # callback is passed the jobId (=Non
      # if not in ccp4i2-db context) and processId that
      # can serve at identifier for subProcess
      import sys
      import time
      # Get the exit status and if successful get the CELL from the outputData
      pid =  ret.get('pid',None)
      #print 'demo_multi_mtzdump.handleDone',ret
      for p in self.subProcessList:
        if p.processId == pid:
          if ret.get('finishStatus') == CPluginScript.SUCCEEDED:
               print(p.container.inputData.HKLIN,p.container.outputData.CELL);sys.stdout.flush()
          self.subProcessList.remove(p)
          break

      if len(self.subProcessList)==0:
        self.reportStatus(CPluginScript.SUCCEEDED)
        print('FINISHED')
      #else:
      #  print len(self.subProcessList),'jobs remaining';sys.stdout.flush()

      return

    @QtCore.Slot()
    def handleTimeout(self):
      #print 'demo_multi_mtzdump.handleTimout', self.subProcessList
      import sys;sys.stdout.flush()
      
      for p in self.subProcessList:
        print('TERMINATING', p.processId,sys.stdout.flush())
        try:
          p.terminate()
        except:
          pass

      self.appendErrorReport(40,str(demo_multi_mtzdump.TIMEOUT_PERIOD))
      self.reportStatus(CPluginScript.FAILED)
