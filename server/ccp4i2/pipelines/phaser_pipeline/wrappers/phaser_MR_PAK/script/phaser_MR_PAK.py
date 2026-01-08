import os
import pickle
import sys

from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script import phaser_MR_AUTO


class phaser_MR_PAK(phaser_MR_AUTO.phaser_MR_AUTO):

    TASKNAME = 'phaser_MR_PAK'                                  # Task name - should be same as class name
    TASKTITLE='Packing function - PHASER'
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    RUNEXTERNALPROCESS=False
    WHATNEXT = ['phaser_expert']

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' },}

    def startProcess(self, command, **kw):
        import phaser
        outputObject = phaser.Output()
        outputObject.setPhenixCallback(self.callbackObject)

        self.prepareCaptureCPlusPlusStdoutToLog()
        resultObject = self.runMR_DAT(outputObject)
        self.finishCaptureCPlusPlusStdout()

        if resultObject == CPluginScript.FAILED: return CPluginScript.FAILED

        self.inputHall = resultObject.getSpaceGroupHall()
        self.inputSpaceGroup = resultObject.getSpaceGroupName()
        inputObject = phaser.InputMR_PAK()
        inputObject.setKILL_FILE(os.path.join(self.getWorkDirectory(),'INTERRUPT'))

        inputObject.setSPAC_HALL(resultObject.getSpaceGroupHall())
        inputObject.setCELL6(resultObject.getUnitCell())
        #print '\n\n\nresutObject',dir(resultObject)
        if self.setKeywords(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseEnsembles(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseSolutions(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.container.inputData.KILLFILEPATH.isSet():
            inputObject.setKILL_FILE(self.container.inputData.KILLFILEPATH.__str__())
        else:
            inputObject.setKILL_FILE(os.path.join(self.getWorkDirectory(),'INTERRUPT'))


        #Now run the main calculation....do something to catch the stdout from the
        #underlying C++ calls
        inputObject.setMUTE(False)
        self.prepareCaptureCPlusPlusStdoutToLog()

        inputObject.setKEYW(True)
        self.resultObject = phaser.runMR_PAK(inputObject, outputObject)

        self.finishCaptureCPlusPlusStdout()
        if not self.resultObject.Success():
            self.appendErrorReport(105, self.resultObject.ErrorName() + '-' + self.resultObject.ErrorMessage())
            return CPluginScript.FAILED
            
        self.analyseResults(self.resultObject)
        return CPluginScript.SUCCEEDED

    # process one or more output files
    # also writes the XML file, previously done by postProcess()
    def processOutputFiles(self):
        import phaser
        resultObject = self.resultObject
        solutions = resultObject.getDotSol()
        if len(solutions) > 0:
            if sys.version_info > (3,0):
                picklePath = str(self.container.outputData.SOLOUT.fullPath)
                with open(picklePath,'wb') as pickleFile:
                    try:
                        pickle.dump(solutions, pickleFile)
                    except:
                        raise
                        print('Unable to Pickle solutions')
                    self.container.outputData.SOLOUT.annotation = 'Solutions from Phaser'
            else:
                picklePath = str(self.container.outputData.SOLOUT.fullPath)
                with open(picklePath,'w') as pickleFile:
                    try:
                        pickle.dump(solutions, pickleFile)
                    except:
                        raise
                        print('Unable to Pickle solutions')
                    self.container.outputData.SOLOUT.annotation = 'Solutions from Phaser'

        #Remove warnings and replace with ones parsed from the resultObject
        if len(self.xmlroot.xpath('PhaserWarnings')) > 0:
            phaser_warnings = [wrng for wrng in resultObject.warnings()]
            for warningsElement in self.xmlroot.xpath('PhaserWarnings')[0]:
                if warningsElement.text not in phaser_warnings:
                    advisoriesElement = etree.SubElement(self.xmlroot,'PhaserAdvisories')
                    advisoryElement = etree.SubElement(advisoriesElement,'Advisory')
                    advisoryElement.text = warningsElement.text
            for warningsElement in self.xmlroot.xpath('PhaserWarnings')[0]:
                warningsElement.getparent().remove(warningsElement)
            for warning in phaser_warnings:
                warningsElement = etree.SubElement(self.xmlroot,'PhaserWarnings')
                warningElement = etree.SubElement(warningsElement,'Warning')
                warningElement.text = warning

        #Remove old digested summaries and add new ones parsed from the result summary block
        for summaryNode in self.xmlroot.xpath('Summary'):
            summaryNode.getparent().remove(summaryNode)
        summary_buffer = '***'
        for text in resultObject.summary().split('\n'):
            if text.startswith("**********") and not summary_buffer.strip().endswith("***"):
                summaryNode = etree.SubElement(self.xmlroot,'Summary')
                summaryNode.text = summary_buffer
                summary_buffer = ""
            summary_buffer += (text+'\n')
        summaryNode = etree.SubElement(self.xmlroot,'Summary')
        summaryNode.text = summary_buffer

        self.flushXML(self.xmlroot)
        return CPluginScript.SUCCEEDED
