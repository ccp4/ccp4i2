import os

from lxml import etree
import phaser

from ......core.CCP4PluginScript import CPluginScript
from ...phaser_MR_AUTO.script import phaser_MR_AUTO


class phaser_MR_FTF(phaser_MR_AUTO.phaser_MR_AUTO):

    TASKNAME = 'phaser_MR_FTF'                                  # Task name - should be same as class name
    TASKTITLE='Translation function - PHASER'
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS=False
    WHATNEXT = ['phaser_expert']

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' },}

    def startProcess(self, command, **kw):
        outputObject = phaser.Output()
        outputObject.setPhenixCallback(self.callbackObject)

        self.prepareCaptureCPlusPlusStdoutToLog()
        resultObject = self.runMR_DAT(outputObject)
        self.finishCaptureCPlusPlusStdout()

        if resultObject == CPluginScript.FAILED: return CPluginScript.FAILED

        self.inputHall = resultObject.getSpaceGroupHall()
        self.inputSpaceGroup = resultObject.getSpaceGroupName()

        inputObject = phaser.InputMR_FTF()

        inputObject.setSPAC_HALL(resultObject.getSpaceGroupHall())
        inputObject.setCELL6(resultObject.getUnitCell())
        if self.container.inputData.F_OR_I.isSet() and self.container.inputData.F_OR_I.__str__() == 'I':
            inputObject.setREFL_I_SIGI(resultObject.getMiller(),resultObject.getIobs(),resultObject.getSigIobs())
        else:
            inputObject.setREFL_F_SIGF(resultObject.getMiller(),resultObject.getFobs(),resultObject.getSigFobs())
        if self.setKeywords(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseContent(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseEnsembles(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.addSearches(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseSolutions(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.container.inputData.KILLFILEPATH.isSet():
            inputObject.setKILL_FILE(self.container.inputData.KILLFILEPATH.__str__())
        else:
            inputObject.setKILL_FILE(os.path.join(self.getWorkDirectory(),'INTERRUPT'))

        inp = self.container.inputData
        if inp.SGALT_SELECT.isSet():
            inputObject.setSGAL_SELE(str(inp.SGALT_SELECT))
            #print 'Setting SGAL_SELE to ',str(inp.SGALT_SELECT)
            if inp.SGALT_SELECT.__str__() == 'LIST' and inp.SGALT_TEST.isSet():
                for sgAltTest in inp.SGALT_TEST:
                    inputObject.addSGAL_TEST(sgAltTest.__str__())

        self.prepareCaptureCPlusPlusStdoutToLog()
        inputObject.setMUTE(False)
        self.resultObject = phaser.runMR_FTF(inputObject, outputObject)
        self.finishCaptureCPlusPlusStdout()

        if not self.resultObject.Success():
            self.appendErrorReport(105, self.resultObject.ErrorName() + '-' + self.resultObject.ErrorMessage())
            return CPluginScript.FAILED
                
        if self.analyseResults(self.resultObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    # process one or more output files
    # also writes the XML file, previously done by postProcess()
    def processOutputFiles(self):
        resultObject = self.resultObject
        solutions = resultObject.getDotSol()
        if len(solutions) > 0:
            picklePath = str(self.container.outputData.SOLOUT.fullPath)
            with open(picklePath,'wb') as pickleFile:
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
    
