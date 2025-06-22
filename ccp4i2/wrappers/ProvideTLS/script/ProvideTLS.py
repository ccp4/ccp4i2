import xml.etree.ElementTree as ET

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class ProvideTLS(CPluginScript):

    TASKNAME = 'ProvideTLS'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        
        self.checkOutputData()

        with open(self.container.outputData.TLSFILE.fullPath.__str__(),"w") as myFile:
            myFile.write(self.container.controlParameters.TLSTEXT.__str__() )

        root = ET.Element('ProvideTLSOutput')
        tlsElement = ET.SubElement(root,'TLSProvided')
        tlsElement.text = self.container.controlParameters.TLSTEXT.__str__()
        CCP4Utils.writeXml(root, self.makeFileName('PROGRAMXML'))

        self.reportStatus(CPluginScript.SUCCEEDED)
