
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Utils


class ProvideTLS(CPluginScript):

    TASKNAME = 'ProvideTLS'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        
        self.checkOutputData()
        
        import shutil
        with open(self.container.outputData.TLSFILE.fullPath.__str__(),"w") as myFile:
            myFile.write(self.container.controlParameters.TLSTEXT.__str__() )
        
        from lxml import etree
        root = etree.Element('ProvideTLSOutput')
        tlsElement = etree.SubElement(root,'TLSProvided')
        tlsElement.text = self.container.controlParameters.TLSTEXT.__str__()
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            CCP4Utils.writeXML(xmlFile,etree.tostring(root,pretty_print=True))

        self.reportStatus(CPluginScript.SUCCEEDED)
