from __future__ import print_function


from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
import pathlib
#from lxml import etree
from xml.etree import ElementTree as ET

class csymmatch(CPluginScript):

    TASKTITLE='Csymmatch - move model to a reference by using symmetry'
    TASKNAME = 'csymmatch'
    TASKMODULE = 'molecular_replacement'
    TASKCOMMAND = 'csymmatch'
    TASKVERSION = 0.0
    ASYNCHRONOUS = False
    MAINTAINER = 'liz.potterton@york.ac.uk'

    def makeCommandAndScript(self):

      inp = self.container.inputData
      par = self.container.controlParameters
      out = self.container.outputData

      import os
      if inp.XYZIN_QUERY.isSelectionSet():
        xyzin_query_file = os.path.join(self.getWorkDirectory(),'XYZIN_QUERY_sel.pdb')
        self.container.inputData.XYZIN_QUERY.loadFile()
        if self.container.inputData.XYZIN_QUERY.isMMCIF():
            xyzin_query_file = str(pathlib.Path(xyzin_query_file).with_suffix('.cif'))
            out.XYZOUT.setFullPath(str(pathlib.Path(out.XYZOUT.fullPath.__str__()).with_suffix('.cif')))
        inp.XYZIN_QUERY.getSelectedAtomsPdbFile(xyzin_query_file)
      else:
        xyzin_query_file = inp.XYZIN_QUERY.fullPath.__str__()
        self.container.inputData.XYZIN_QUERY.loadFile()
        if self.container.inputData.XYZIN_QUERY.isMMCIF():
            out.XYZOUT.setFullPath(str(pathlib.Path(out.XYZOUT.fullPath.__str__()).with_suffix('.cif')))
      if inp.XYZIN_TARGET.isSelectionSet():
        xyzin_target_file = os.path.join(self.getWorkDirectory(),'XYZIN_TARGET_sel.pdb')
        self.container.inputData.XYZIN_TARGET.loadFile()
        if self.container.inputData.XYZIN_TARGET.isMMCIF():
            xyzin_target_file = str(pathlib.Path(xyzin_target_file).with_suffix('.cif'))
        inp.XYZIN_TARGET.getSelectedAtomsPdbFile(xyzin_target_file)
      else:
        xyzin_target_file = inp.XYZIN_TARGET.fullPath.__str__()
      self.appendCommandLine( [ '-pdbin', xyzin_query_file ] )
      self.appendCommandLine( [ '-pdbin-ref', xyzin_target_file ] )
      if par.ORIGIN_HAND.isSet():
        if par.ORIGIN_HAND:
          self.appendCommandLine( [ '-origin-hand'] )
      if par.CONNECTIVITY_RADIUS.isSet():
        self.appendCommandLine( [ '-connectivity-radius',str(par.CONNECTIVITY_RADIUS)] )
      self.appendCommandLine( [ '-pdbout' , out.XYZOUT.fullPath.__str__() ] )

    def processOutputFiles(self):
        logName = self.makeFileName('LOG')
        
        xmlRoot = ET.Element('Csymmatch')
        segmentNode = None
        with open (logName,'r') as logFile:
            lines = logFile.readlines()
            for line in lines:
                if line.strip().startswith('Change of hand'):
                    handNode = ET.SubElement(xmlRoot,'ChangeOfHand')
                    handNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('Change of origin'):
                    originNode = ET.SubElement(xmlRoot,'ChangeOfOrigin')
                    originNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('Chain'):
                    segmentNode = ET.SubElement(xmlRoot,'Segment')
                    rangeNode = ET.SubElement(segmentNode,'Range')
                    rangeNode.text = line.strip().split('will')[0]
                elif line.strip().startswith('Symmetry operator'):
                    if segmentNode is not None:
                        operatorNode = ET.SubElement(segmentNode,'Operator')
                        operatorNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('Lattice shift'):
                    if segmentNode is not None:
                        shiftNode = ET.SubElement(segmentNode,'Shift')
                        shiftNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('with normalised score'):
                    if segmentNode is not None:
                        scoreNode = ET.SubElement(segmentNode,'Score')
                        scoreNode.text = line.strip().split(':')[1]
    
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            xmlString = ET.tostring(xmlRoot)
            CCP4Utils.writeXML(xmlFile,xmlString)


        return CPluginScript.SUCCEEDED
#---------------------------------------------
import unittest

class testcsymmatch( unittest.TestCase ) :

   def setUp(self):
    from core import CCP4Modules
    self.app = CCP4Modules.QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from core import CCP4Modules
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)

   def test1( self ) :

      from core import CCP4Modules
      import os

      workDirectory = CCP4Utils.getTestTmpDir()
      xmlInput = os.path.join( CCP4Utils.getCCP4I2Dir(), 'wrappers', 'csymmatch', 'test_data', 'test1'+'.params.xml' )
      self.wrapper = csymmatch(parent=CCP4Modules.QTAPPLICATION(), name='csymmatch_test1',workDirectory=workDirectory)
      self.wrapper.container.loadDataFromXml( xmlInput )

      self.wrapper.setWaitForFinished( 1000000 )
      pid = self.wrapper.process()
      self.wrapper.setWaitForFinished( -1 )
      if len(self.wrapper.errorReport)>0:
         print(self.wrapper.errorReport.report())

def TESTSUITE() :

   suite = unittest.TestLoader().loadTestsFromTestCase( testcsymmatch )
   return suite

def testModule() :

   suite = TESTSUITE()
   unittest.TextTestRunner( verbosity=2 ).run( suite )

