from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript


class pisapipe(CPluginScript):

    TASKMODULE = 'test'                              # Where this plugin will appear on the gui
    TASKTITLE = ' Structure analysis with Pisa'     # A short title for gui menu
    TASKNAME = 'pisapipe'                                  # Task name - should be same as class name

    def process(self):
        identifier = str(self.container.inputData.PDBIN.dbFileId)
        btask = self.makePluginObject('pisa_analyse')
        btask.container.controlParameters.IDENTIFIER = identifier
        btask.container.inputData.PDBIN = self.container.inputData.PDBIN
        jobStatus = btask.process()
        if jobStatus == CPluginScript.FAILED:
            self.reportStatus(jobStatus)
            return

        self.xmlTask = self.makePluginObject('pisa_xml')
        self.xmlTask.container.controlParameters.IDENTIFIER = identifier
        status = self.xmlTask.process()
        print('pisa_xmlFinished', status)
        if status == CPluginScript.FAILED:
            self.reportStatus(status)

        from ccp4i2.core import CCP4Utils
        self.xmlroot = etree.Element('pisapipe')
        xmlOfTask = CCP4Utils.openFileToEtree(self.xmlTask.makeFileName('PROGRAMXML'))
        self.xmlroot.append(xmlOfTask)
        with open(self.makeFileName('PROGRAMXML'),'w') as outputXMLFile:
            CCP4Utils.writeXML(outputXMLFile, etree.tostring(self.xmlroot,pretty_print=True))
        self.reportStatus(CPluginScript.SUCCEEDED)
