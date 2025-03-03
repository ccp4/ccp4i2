from __future__ import print_function


from core import CCP4PluginScript
from core import CCP4XtalData

  
class chltofom(CCP4PluginScript.CPluginScript):


    TASKNAME = 'chltofom'                                  # Task name - should be same as class name
    TASKCOMMAND = 'chltofom'                             # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    MAINTAINER = 'liz.potterton@york.ac.uk'
    PURGESEARCHLIST = [  [ 'hklout.mtz' , 0 ],    
                         [ 'HKLOUT.mtz' , 0 ]
                       ]
    
    def makeCommandAndScript(self):
      self.appendCommandLine( ['-mtzin',self.container.inputData.HKLIN.__str__() ] )
      if self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM:
        self.container.controlParameters.DIRECTION.set('FOMTOHL') 
        self.appendCommandLine ( ['-colin-phifom', '/*/*/[PHI,FOM]' ])
        self.appendCommandLine ( ['-colout', '/*/*/[HLA,HLB,HLC,HLD]' ])
      else:
        self.appendCommandLine ( ['-colin-hl', '/*/*/[HLA,HLB,HLC,HLD]' ])
        self.appendCommandLine ( ['-colout','/*/*/[PHI,FOM]' ])
      if self.container.controlParameters.OUTPUTMINIMTZ:
          from core import CCP4Utils
          self.tmpHklout = CCP4Utils.makeTmpFile(extension='mtz')
          self.appendCommandLine( [ '-mtzout',self.tmpHklout] )
      else:
        self.appendCommandLine( [ '-mtzout',self.container.outputData.HKLOUT.__str__()] )
      return CCP4PluginScript.CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      #print 'chltofom.processOutputFiles',self.container.controlParameters.OUTPUTMINIMTZ,self.__dict__.get('tmpHklout','NONE')
      if self.container.controlParameters.OUTPUTMINIMTZ:
        import os
        logFile = os.path.splitext(self.tmpHklout)[0]+'.log'

        if self.container.inputData.HKLIN.annotation.isSet():
          anno = str(self.container.inputData.HKLIN.annotation)
        else:
          anno = self.container.outputData.HKLOUT.qualifiers('guiLabel')
       
        if self.container.controlParameters.DIRECTION.isSet() and self.container.controlParameters.DIRECTION.__str__() == 'FOMTOHL':
          status = self.splitMtz(self.tmpHklout,[[self.container.outputData.HKLOUT.__str__(),'HLA,HLB,HLC,HLD']],logFile)
          self.container.outputData.HKLOUT.annotation = anno + ' as HLcoeffs'
        else:
          status = self.splitMtz(self.tmpHklout,[[self.container.outputData.HKLOUT.__str__(),'PHI,FOM']],logFile)
          self.container.outputData.HKLOUT.annotation = anno  + ' as phi/FOM'
          self.container.outputData.HKLOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM
        if status != CCP4PluginScript.CPluginScript.SUCCEEDED: return status
        if not os.path.exists(self.container.outputData.HKLOUT.__str__()):
          try:
            shutil.copyfile(self.tmpHklfile,self.container.outputData.HKLOUT.__str__())
          except:
            pass
          return CCP4PluginScript.CPluginScript.FAILED
      return CCP4PluginScript.CPluginScript.SUCCEEDED
     
#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testchltofom(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
    self.app = QTAPPLICATION()

   def tearDown(self):
    from core.CCP4Modules import PROCESSMANAGER
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     from core.CCP4Modules import QTAPPLICATION
     from core import CCP4Utils
     import os

     workDirectory = CCP4Utils.getTestTmpDir()
     logFile = os.path.join(workDirectory,'chltofom_test1.log')
     # Delete any existing log file
     if os.path.exists(logFile): os.remove(logFile)
     outFile = os.path.join(workDirectory,'chltofom_test1.mtz')
     print('testchltofom outFile',outFile)
     if os.path.exists(outFile): os.remove(outFile)
     wrapper = chltofom(parent=QTAPPLICATION(),name='chltofom_test1')
     wrapper.container.inputData.HKLIN.setFullPath(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','rnase25_mini_HL.mtz'))
     wrapper.container.outputData.HKLOUT.setFullPath(outFile)
     wrapper.container.controlParameters.OUTPUTMINIMTZ.set(True)
     pid = wrapper.process()
     self.assertTrue(os.path.exists(  outFile),'No output file from chltofom_test1')                                     
     wrapper.container.outputData.HKLOUT.loadFile()
     columns = wrapper.container.outputData.HKLOUT.fileContent.getListOfColumns()
     print('chltofom.processOutputFiles',columns)
     self.assertEqual(len(columns),2,'Output from chltofom_test1 has wrong number of columns')
     self.assertTrue(columns[0].columnLabel.__str__() in ['PHI','FOM'] and columns[1].columnLabel.__str__() in ['PHI','FOM'],'Output from chltofom_test1 has wrong column labels')


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testchltofom)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
