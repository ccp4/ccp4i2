
import os
import xml.etree.ElementTree as ET

from PySide2 import QtCore

from ......core.CCP4PluginScript import CPluginScript
from ......core import CCP4Utils
from ......core.CCP4Modules import LAUNCHER


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
        self.xmlroot = ET.Element('pisa_xml')
        result = self.retrieveAssemblies()
        if result != CPluginScript.SUCCEEDED:
            self.reportStatus(result)

    def retrieveAssemblies(self):
        try:
            self.assembliesXMLPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'assemblies.xml'))
            LAUNCHER().launch(viewer='pisa', argList = [self.container.controlParameters.IDENTIFIER.__str__(), '-xml','assemblies'], callBack = self.retrieveAssembliesFinished, logFile = self.assembliesXMLPath)
        except:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    @QtCore.Slot(int,int)
    def retrieveAssembliesFinished(self,exitCode,exitStatus):
        assembliesNode = ET.SubElement(self.xmlroot,'Assemblies')
        assembliesXML = ET.parse(self.assembliesXMLPath).getroot()
        assembliesNode.append(assembliesXML)
        result = self.retrieveInterfaces()
        if result != CPluginScript.SUCCEEDED:
            self.reportStatus(result)

    def retrieveInterfaces(self):
        try:
            self.interfacesXMLPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'interfaces.xml'))
            LAUNCHER().launch(viewer='pisa', argList = [self.container.controlParameters.IDENTIFIER.__str__(), '-xml','interfaces'], callBack = self.retrieveInterfacesFinished, logFile = self.interfacesXMLPath)
        except:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    @QtCore.Slot(int,int)
    def retrieveInterfacesFinished(self,exitCode,exitStatus):
        interfacesNode = ET.SubElement(self.xmlroot,'Interfaces')
        interfacesXML = ET.parse(self.interfacesXMLPath).getroot()
        interfacesNode.append(interfacesXML)
        self.finishWithStatus(CPluginScript.SUCCEEDED)

    def flushXML(self):
        CCP4Utils.writeXml(self.xmlroot, self.makeFileName('PROGRAMXML'))

    def finishWithStatus(self, status=CPluginScript.SUCCEEDED):
        self.flushXML()
        self.reportStatus(status)
