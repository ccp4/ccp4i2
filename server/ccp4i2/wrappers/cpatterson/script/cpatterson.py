from lxml import etree

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class cpatterson(CPluginScript):

    TASKMODULE = 'wrappers'
    TASKTITLE = 'Prepare map coefficients'
    TASKNAME = 'cpatterson'
    TASKCOMMAND = 'cpatterson'
    TASKVERSION= 0.0
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    def processInputFiles ( self ):
        from ccp4i2.core import CCP4XtalData
        self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])

    def makeCommandAndScript(self):

        out = self.container.outputData

        self.appendCommandLine([ '-stdin' ])

        self.appendCommandScript([ '-mtzin ' + self.hklin ])
        self.appendCommandScript([ '-mapout %s'%(str(out.MAPOUT.fullPath)) ])
        self.appendCommandScript ([ '-colin-fo F,SIGF' ])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        self.container.outputData.MAPOUT.annotation = 'Computed using ' + str(self.container.inputData.F_SIGF.annotation)

        xmlRoot = etree.Element('Patterson')
        with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
            xmlString = etree.tostring ( xmlRoot, pretty_print=True )
            CCP4Utils.writeXML(xmlFile,xmlString)

        return CPluginScript.SUCCEEDED
