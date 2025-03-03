from __future__ import print_function



from core.CCP4PluginScript import CPluginScript
from wrappers.x2mtz.script import x2mtz

class scalepack2mtz(x2mtz.x2mtz):

    TASKMODULE = 'test'
    TASKTITLE = 'Convert scalepack merged reflection file to MTZ'
    TASKNAME = 'scalepack2mtz'  
    TASKVERSION= 0.1    
    TASKCOMMAND = 'scalepack2mtz'  
    MAINTAINER = 'liz.potterton@york.ac.uk'

    def makeCommandAndScript(self):

      inp = self.container.inputData
      # input parameters
      self.appendCommandLine(["HKLIN", inp.HKLIN.fullPath])

      if self.container.inputData.SPACEGROUPCELL.cell.isSet():
         inp.SPACEGROUPCELL.cell.set(inp.SPACEGROUPCELL.cell.fix(inp.SPACEGROUPCELL.cell.get()))
         self.appendCommandScript("CELL %f %f %f %f %f %f" %
           (inp.SPACEGROUPCELL.cell.a.__float__(), inp.SPACEGROUPCELL.cell.b.__float__(), inp.SPACEGROUPCELL.cell.c.__float__(),
            inp.SPACEGROUPCELL.cell.alpha.__float__(), inp.SPACEGROUPCELL.cell.beta.__float__(), inp.SPACEGROUPCELL.cell.gamma.__float__()))
      if inp.SPACEGROUPCELL.spaceGroup.isSet():
         self.appendCommandScript("SYMMETRY %d" % inp.SPACEGROUPCELL.spaceGroup.number())     
      if inp.WAVELENGTH.isSet():
          self.appendCommandScript("WAVE %f" % inp.WAVELENGTH.__float__() )                                                
      if inp.CRYSTALNAME.isSet():self.appendCommandScript("XNAME '%s'"%str(inp.CRYSTALNAME))
      if inp.DATASETNAME.isSet():self.appendCommandScript("DNAME '%s'"%str(inp.DATASETNAME))

      # control parameters
      if self.container.controlParameters.ANOMALOUS:
          self.appendCommandScript("anomalous")

      s = self.resolutionRangeCommand()
      if s != "":
          self.appendCommandScript(s)
          
      # output parameters
      if self.container.outputData.HKLOUT.isSet():
         self.appendCommandLine(["HKLOUT", self.container.outputData.HKLOUT.fullPath])
      self.appendCommandScript('END')

      return 0

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def resolutionRangeCommand(self):
        r1 = self.container.inputData.RESOLUTION_RANGE.start
        r2 = self.container.inputData.RESOLUTION_RANGE.end
        print("S2M resolutionRangeCommand",
              self.container.inputData.RESOLUTION_RANGE,
              r1, r2)
        
        high = 0.0
        low  = 99999.0
        if not r1.isSet() and not r2.isSet():
            # nothing set
            s = ""
        else:
            if r1.isSet():
                low = float(r1)
            if r2.isSet():
                high = float(r2)
            if low < high:
                low, high = high, low
            s = "RESOLUTION %f  %f" % (low, high)

        return s


#======================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testscalepack2mtz(unittest.TestCase):

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
     from core import CCP4Modules
     from core import CCP4Utils
     import os

     workDirectory = CCP4Utils.getTestTmpDir()
     # this needs to agree with name attribute below
     logFile = os.path.join(workDirectory,'convert2mtz_test_auto.log')
     # Delete any existing log file
     if os.path.exists(logFile): os.remove(logFile)

     self.wrapper = scalepack2mtz(parent=CCP4Modules.QTAPPLICATION(),name='convert2mtz_test_auto',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','convert2mtz','test_data','cns_input.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
     #self.assertTrue(os.path.exists(logFile),'No log file found')


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testconvert2mtz)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
