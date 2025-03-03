from __future__ import print_function

"""
     chainsaw.py: CCP4 GUI Project
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

import os
from core.CCP4PluginScript import CPluginScript

class chainsaw(CPluginScript):

    TASKMODULE = 'molecular replacement' # Where this plugin will appear on gui
    TASKTITLE = 'edit search model' # A short title for gui menu
    TASKNAME = 'chainsaw'   # Task name - should be same as class name
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'ronan.keegan@stfc.ac.uk'
    PERFORMANCECLASS = 'CAtomCountPerformance'

    # used by the base class startProcess()
    TASKCOMMAND = 'chainsaw'   # The command to run the executable

    def makeCommandAndScript(self):
      print('chainsaw.makeCommandAndScript')

      if self.container.inputData.XYZIN.isSelectionSet():
        xyzin_file = os.path.join(self.getWorkDirectory(),"XYZIN_sel.pdb")
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(xyzin_file)
      else:
        xyzin_file = str(self.container.inputData.XYZIN)
      self.appendCommandLine(['XYZIN',xyzin_file])
      self.appendCommandLine(['ALIGNIN',self.inputAlignmentFileName])
      self.appendCommandLine(['XYZOUT',self.container.outputData.XYZOUT.fullPath])

      ### keywords
      self.appendCommandScript('MODE %s' % self.container.controlParameters.MODE)
      self.appendCommandScript('END')

    def processInputFiles(self):
        self.cryst1card = None
        with open(str(self.container.inputData.XYZIN.fullPath),'r') as inputFile:
            lines = inputFile.readlines()
            for line in lines:
                if line.startswith('CRYST1'):
                    self.cryst1card = line

        # Ensure correct alignment file extension and reorder file if target is not first
        formt,idList = self.container.inputData.ALIGNIN.identifyFile()
        ext = { 'unknown' : 'aln',
                'clustal' : 'aln' ,
                'fasta' : 'fas',
                'pir' : 'pir' ,
                'phylip' : 'phy' }.get(formt,'aln')
        self.inputAlignmentFileName = os.path.normpath(os.path.join(self.getWorkDirectory(),'tempAlignment.'+ext))
        #print 'inputAlignmentFileName',self.inputAlignmentFileName,'TARGETINDEX',self.container.controlParameters.TARGETINDEX
        if self.container.controlParameters.TARGETINDEX.isSet() and self.container.controlParameters.TARGETINDEX != 0:
          formt,idList = self.container.inputData.ALIGNIN.identifyFile()
          self.container.inputData.ALIGNIN.convertFormat(formt,self.inputAlignmentFileName,reorder='reverse')
        else:
          # Forcing the file extension to chainsaw requirement
          import shutil
          shutil.copyfile(self.container.inputData.ALIGNIN.__str__(), self.inputAlignmentFileName)
          return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      if self.cryst1card is not None:
          import os
          tmpFilename = str(self.container.outputData.XYZOUT.fullPath)+'_tmp'
          os.rename(str(self.container.outputData.XYZOUT.fullPath), tmpFilename)
          with open(tmpFilename,'r') as inputFile:
              lines = inputFile.readlines()
              with open(str(self.container.outputData.XYZOUT.fullPath),'w') as outputFile:
                  outputFile.write(self.cryst1card)
                  for line in lines:
                      outputFile.write(line)
      self.container.outputData.XYZOUT.subType=2

      self.container.outputData.XYZOUT.annotation = 'Chainsawed model coordinates'

      # sanity check that chainsaw has produced something
      # -1 means we have not managed to get value out of log file
      self.container.outputData.NO_DELETED = -1
      self.container.outputData.NO_CONSERVED = -1
      self.container.outputData.NO_MUTATED = -1
      logText = self.logFileText()
      pyListLogLines = logText.split("\n")
      for j, pyStrLine in enumerate(pyListLogLines):
         if "conserved" in pyStrLine and "mutated" in pyStrLine:
            self.container.outputData.NO_DELETED = int(pyStrLine.split()[0])
            self.container.outputData.NO_CONSERVED = int(pyStrLine.split()[2])
            self.container.outputData.NO_MUTATED = int(pyStrLine.split()[4])

      if self.container.outputData.NO_CONSERVED > -1 and self.container.outputData.NO_MUTATED > -1:
         self.container.outputData.NUMRES_MODEL = self.container.outputData.NO_CONSERVED + self.container.outputData.NO_MUTATED
         self.container.outputData.SEQID = float(self.container.outputData.NO_CONSERVED) / float(self.container.outputData.NUMRES_MODEL)

      self.container.saveDataToXml(self.makeFileName( 'PROGRAMXML' ))

      # Set a performance parameter mostly used in i2 testing
      self.container.outputData.PERFORMANCE.setFromPdbDataFile(self.container.outputData.XYZOUT)
      #print 'PERFORMANCE',self.container.outputData.PERFORMANCE.nAtoms,self.container.outputData.PERFORMANCE.nResidues

      return CPluginScript.SUCCEEDED
    
#======================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testchainsaw(unittest.TestCase):

   def setUp(self):
    from core import CCP4Modules
    self.app = CCP4Modules.QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from core import CCP4Modules
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     from core import CCP4Modules, CCP4Utils
     import os

     workDirectory = CCP4Utils.getTestTmpDir()
     # this needs to agree with name attribute below
     logFile = os.path.join(workDirectory,'chainsaw_test1.log')
     # Delete any existing log file
     if os.path.exists(logFile): os.remove(logFile)

     self.wrapper = chainsaw(parent=CCP4Modules.QTAPPLICATION(),name='chainsaw_test1',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','chainsaw','test_data','chainsaw_test1.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
     #self.assertTrue(os.path.exists(logFile),'No log file found')
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testchainsaw)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
