"""
Copyright (C) 2010 University of York
"""

import xml.etree.ElementTree as ET

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


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

        CCP4Utils.writeXml(xmlRoot, self.makeFileName('PROGRAMXML'))

        return CPluginScript.SUCCEEDED
