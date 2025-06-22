import os
import xml.etree.ElementTree as ET

from lxml import etree

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class AcedrgLink(CPluginScript):
    TASKNAME = 'AcedrgLink'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'nicholls@mrc-lmb.cam.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="acedrg"
    
    def __init__(self, *args, **kws):
        super(AcedrgLink, self).__init__(*args, **kws)

    def processInputFiles(self):
        print('AceDRG in link mode - processing input')
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        self.appendCommandLine(['-L',self.container.inputData.INSTRUCTION_FILE.fullPath])
        self.appendCommandLine(['-o',os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__())])
        if self.container.controlParameters.EXTRA_ACEDRG_KEYWORDS.isSet():
           for kwLine in str(self.container.controlParameters.EXTRA_ACEDRG_KEYWORDS).split('\n'):
              kw = kwLine.lstrip().rstrip()
              #print 'kw','['+str(kw)+']'
              if len(kw)>0:
                 if str(kw)[0] != '#':
                    self.appendCommandLine(kw)
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        print('AceDRG in link mode - processing output')

        self.container.outputData.CIF_OUT.setFullPath(os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__()+"_link.cif"))
        if self.container.inputData.ANNOTATION.isSet():
           self.container.outputData.CIF_OUT.annotation.set('Link dictionary: '+self.container.inputData.ANNOTATION.__str__())
        else:
           self.container.outputData.CIF_OUT.annotation.set('Link dictionary: '+self.container.inputData.LINK_ID.__str__())

        self.container.outputData.UNL_PDB.setFullPath(os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__()+"_TMP","UNL_for_link.pdb"))
        self.container.outputData.UNL_CIF.setFullPath(os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__()+"_TMP","UNL_for_link.cif"))

        #Create (dummy) PROGRAMXML, which basically contains only the log text of the job
        #without this, a report will not be generated
        
        xmlStructure = ET.Element("acedrg_link")
        logText = ET.SubElement(xmlStructure,"LogText")
        with open(self.makeFileName("LOG"),"r") as logFile:
            logText.text = etree.CDATA(logFile.read())
        CCP4Utils.writeXml(xmlStructure, self.makeFileName("PROGRAMXML"))

        return CPluginScript.SUCCEEDED
