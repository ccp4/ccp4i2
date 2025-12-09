
from core.CCP4PluginScript import CPluginScript
from lxml import etree
import os
from core import CCP4Modules
from core import CCP4Utils

from ccp4i2.baselayer import QtCore

class pisa_xml(CPluginScript):

    TASKMODULE = None                               # Where this plugin will appear on the gui
    TASKTITLE = 'Structure analysis with Pisa'     # A short title for gui menu
    TASKNAME = 'pisa_xml'                                  # Task name - should be same as class name
    TASKCOMMAND = 'pisa'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    ASYNCHRONOUS=True

  
    def process(self):
        self.xmlroot = etree.Element('pisa_xml')
        result = self.retrieveAssemblies()
        if result != CPluginScript.SUCCEEDED:
            self.reportStatus(result)

    def retrieveAssemblies(self):
        try:
            self.assembliesXMLPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'assemblies.xml'))
            CCP4Modules.LAUNCHER().launch(viewer='pisa', argList = [self.container.controlParameters.IDENTIFIER.__str__(), '-xml','assemblies'], callBack = self.retrieveAssembliesFinished, logFile = self.assembliesXMLPath)
        except:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    @QtCore.Slot(int,int)
    def retrieveAssembliesFinished(self,exitCode,exitStatus):
        assembliesNode = etree.SubElement(self.xmlroot,'Assemblies')
        assembliesXML = CCP4Utils.openFileToEtree(self.assembliesXMLPath)
        assembliesNode.append(assembliesXML)
        result = self.retrieveInterfaces()
        if result != CPluginScript.SUCCEEDED:
            self.reportStatus(result)

    def retrieveInterfaces(self):
        try:
            self.interfacesXMLPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'interfaces.xml'))
            CCP4Modules.LAUNCHER().launch(viewer='pisa', argList = [self.container.controlParameters.IDENTIFIER.__str__(), '-xml','interfaces'], callBack = self.retrieveInterfacesFinished, logFile = self.interfacesXMLPath)
        except:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    @QtCore.Slot(int,int)
    def retrieveInterfacesFinished(self,exitCode,exitStatus):
        interfacesNode = etree.SubElement(self.xmlroot,'Interfaces')
        interfacesXML = CCP4Utils.openFileToEtree(self.interfacesXMLPath)
        interfacesNode.append(interfacesXML)
        self.finishWithStatus(CPluginScript.SUCCEEDED)

    def flushXML(self):
        with open(self.makeFileName('PROGRAMXML'),'w') as outputXML:
            CCP4Utils.writeXML(outputXML,etree.tostring(self.xmlroot,pretty_print=True))

    def finishWithStatus(self, status=CPluginScript.SUCCEEDED):
        self.flushXML()
        self.reportStatus(status)
