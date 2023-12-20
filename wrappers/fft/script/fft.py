"""
     fft.scripts.py: CCP4 GUI Project
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

from core.CCP4PluginScript import CPluginScript
#from lxml import etree
from xml.etree import ElementTree as ET
from core import CCP4Utils

class fft(CPluginScript):

    TASKMODULE = 'test'
    TASKTITLE = 'Export map'
    TASKNAME = 'fft'
    TASKCOMMAND = 'cfft'
    TASKVERSION= 0.0

    def makeCommandAndScript(self):
        inp = self.container.inputData
        out =  self.container.outputData
        con = self.container.controlParameters

        self.appendCommandLine([ '-stdin' ])

        self.appendCommandScript([ '-mtzin %s'%(str(inp.FPHIIN.fullPath)) ])
        self.appendCommandScript([ '-mapout %s'%(str(out.MAPOUT.fullPath)) ])
        
        if inp.INDEX_U.isSet():
            if inp.INDEX_V.isSet():
                if inp.INDEX_W.isSet():
                    self.appendCommandScript ([ '-grid %s,%s,%s'%(str(inp.INDEX_U))%(str(inp.INDEX_V))%(str(inp.INDEX_W)) ])
        
        if inp.RESOLUTION.isSet():
            self.appendCommandScript ([ '-resolution %s'%(str(inp.RESOLUTION)) ])

        if con.B_VALUE.isSet():
            self.appendCommandScript ([ '-b-value %s'%(str(con.B_VALUE)) ])

        if con.U_VALUE.isSet():
            self.appendCommandScript ([ '-u-value %s'%(str(con.U_VALUE)) ])

        self.appendCommandScript([ '-colin-fc F,PHI' ])
        self.appendCommandScript([ '-stats ' ])    

        return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):
        self.container.outputData.MAPOUT.annotation = 'Computed using ' + str(self.container.inputData.FPHIIN.annotation)
      
        lines = open(self.makeFileName('LOG')).readlines()
        xmlRoot = ET.Element('Cfft')

        readingStuff = False

        for line in lines :
            if len( line.strip().split ( ) ) < 2 :
                readingStuff = False

            if readingStuff :
                if line.strip().startswith ( 'Number of' ) :
                    npoints = ET.SubElement(xmlRoot,'NPoints')
                    npoints.text = line.strip().split ( )[5]
                elif line.strip().startswith ( '1st' ) :
                    firstMomZero = ET.SubElement(xmlRoot, 'FirstMomZero' )
                    firstMomZero.text = line.strip().split ( )[5]
                    firstMomMean = ET.SubElement(xmlRoot, 'FirstMomMean' )
                    firstMomMean.text = line.strip().split ( )[9] 
                elif line.strip().startswith ( '2nd' ) :
                    secondMomZero = ET.SubElement(xmlRoot, 'SecondMomZero' )
                    secondMomZero.text = line.strip().split ( )[5]
                    secondMomMean = ET.SubElement(xmlRoot, 'SecondMomMean' )
                    secondMomMean.text = line.strip().split ( )[9] 
                elif line.strip().startswith ( '3rd' ) :
                    thirdMomZero = ET.SubElement(xmlRoot, 'ThirdMomZero' )
                    thirdMomZero.text = line.strip().split ( )[5]
                    thirdMomMean = ET.SubElement(xmlRoot, 'ThirdMomMean' )
                    thirdMomMean.text = line.strip().split ( )[9] 
                elif line.strip().startswith ( 'Range' ) :
                    minimum = ET.SubElement(xmlRoot, 'Min' )
                    minimum.text = line.strip().split ( )[2]
                    maximum = ET.SubElement(xmlRoot, 'Max' )
                    maximum.text = line.strip().split ( )[4] 
            elif line.strip().startswith ( 'Map statistics' ) :
                readingStuff = True

        with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
            xmlString = ET.tostring ( xmlRoot )
            CCP4Utils.writeXML(xmlFile,xmlString)

        return CPluginScript.SUCCEEDED

#=====================================================================================================
#=================================test suite=========================================================
#=====================================================================================================

import unittest
from core.CCP4Utils import getCCP4I2Dir,getTMP

# unit testing asynchronous processes potential tricky but QProcess has option to wait for finished
 
class testFft(unittest.TestCase):
  
  def setUp(self):
    # make all background jobs wait for completion
    PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)


  def testFft(self):
    import os
    inputData =  CScriptDataContainer(name='fft_test',containerType='inputData',initialise=fft.INPUTDATA)
    outputData =  CScriptDataContainer(name='fft_test',containerType='outputData',initialise=fft.OUTPUTDATA)
    try:
      inputData.importXML(os.path.join(getCCP4I2Dir(),'wrappers','fft','test_data','fft_test_1.def.xml'))
    except CException as e:
      self.fail(e.errorType)
    try:
      outputData.importXML(os.path.join(getCCP4I2Dir(),'wrappers','fft','test_data','fft_test_1.def.xml'))
    except CException as e:
      self.fail(e.errorType)
      
    wrapper = fft()
    pid = wrapper.process()


def testSuite():
  suite = unittest.TestLoader().loadTestsFromTestCase(testFft)
  return suite

def runAllTests():
  suite = testSuite()
  unittest.TextTestRunner(verbosity=2).run(suite)
