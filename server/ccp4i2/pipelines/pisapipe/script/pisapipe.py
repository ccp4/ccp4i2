from ccp4i2.baselayer import QtCore
from ccp4i2.core.CCP4PluginScript import CPluginScript


class pisapipe(CPluginScript):

    TASKMODULE = 'test'                              # Where this plugin will appear on the gui
    TASKTITLE = ' Structure analysis with Pisa'     # A short title for gui menu
    TASKNAME = 'pisapipe'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template


    def process(self):
      #identifier = str(self.container.paramsHeader.project) +'_'+str(self.container.paramsHeader.jobNumber)
      identifier = str(self.container.inputData.PDBIN.dbFileId)
      btask = self.makePluginObject('pisa_analyse')
      btask.container.controlParameters.IDENTIFIER = identifier
      # set synchronous job control
      btask.setWaitForFinished(1000000)
      btask.container.inputData.PDBIN = self.container.inputData.PDBIN
      jobStatus = btask.process()
      if jobStatus == CPluginScript.FAILED:
          self.reportStatus(jobStatus)
          return
      self.xmlTask = self.makePluginObject('pisa_xml')
      self.xmlTask.container.controlParameters.IDENTIFIER = identifier
      self.connectSignal(self.xmlTask,'finished',self.pisa_xmlFinished)
      self.xmlTask.process()

    @QtCore.Slot(dict)
    def pisa_xmlFinished(self, statusDict):
        status = statusDict['finishStatus']
        print('pisa_xmlFinished', status)
        if status == CPluginScript.FAILED: self.reportStatus(status)
        from lxml import etree

        from ccp4i2.core import CCP4Utils
        self.xmlroot = etree.Element('pisapipe')
        xmlOfTask = CCP4Utils.openFileToEtree(self.xmlTask.makeFileName('PROGRAMXML'))
        self.xmlroot.append(xmlOfTask)
        with open(self.makeFileName('PROGRAMXML'),'w') as outputXMLFile:
            CCP4Utils.writeXML(outputXML,etree.tostring(self.xmlroot,pretty_print=True))
        self.reportStatus(CPluginScript.SUCCEEDED)
