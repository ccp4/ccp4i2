
import os

from lxml import etree

from ccp4i2.baselayer import QtCore
from ccp4i2.core import CCP4Modules, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class pisa_xml(CPluginScript):

    TASKTITLE = 'Structure analysis with Pisa'
    TASKNAME = 'pisa_xml'
    TASKCOMMAND = 'pisa'

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
