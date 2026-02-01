from lxml import etree

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class fft(CPluginScript):

    TASKMODULE = 'test'
    TASKTITLE = 'Export map'
    TASKNAME = 'fft'
    TASKCOMMAND = 'cfft'

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
        xmlRoot = etree.Element('Cfft')

        readingStuff = False

        for line in lines :
            if len( line.strip().split ( ) ) < 2 :
                readingStuff = False

            if readingStuff :
                if line.strip().startswith ( 'Number of' ) :
                    npoints = etree.SubElement(xmlRoot,'NPoints')
                    npoints.text = line.strip().split ( )[5]
                elif line.strip().startswith ( '1st' ) :
                    firstMomZero = etree.SubElement(xmlRoot, 'FirstMomZero' )
                    firstMomZero.text = line.strip().split ( )[5]
                    firstMomMean = etree.SubElement(xmlRoot, 'FirstMomMean' )
                    firstMomMean.text = line.strip().split ( )[9] 
                elif line.strip().startswith ( '2nd' ) :
                    secondMomZero = etree.SubElement(xmlRoot, 'SecondMomZero' )
                    secondMomZero.text = line.strip().split ( )[5]
                    secondMomMean = etree.SubElement(xmlRoot, 'SecondMomMean' )
                    secondMomMean.text = line.strip().split ( )[9] 
                elif line.strip().startswith ( '3rd' ) :
                    thirdMomZero = etree.SubElement(xmlRoot, 'ThirdMomZero' )
                    thirdMomZero.text = line.strip().split ( )[5]
                    thirdMomMean = etree.SubElement(xmlRoot, 'ThirdMomMean' )
                    thirdMomMean.text = line.strip().split ( )[9] 
                elif line.strip().startswith ( 'Range' ) :
                    minimum = etree.SubElement(xmlRoot, 'Min' )
                    minimum.text = line.strip().split ( )[2]
                    maximum = etree.SubElement(xmlRoot, 'Max' )
                    maximum.text = line.strip().split ( )[4] 
            elif line.strip().startswith ( 'Map statistics' ) :
                readingStuff = True

        with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
            xmlString = etree.tostring ( xmlRoot, pretty_print=True )
            CCP4Utils.writeXML(xmlFile,xmlString)

        return CPluginScript.SUCCEEDED
