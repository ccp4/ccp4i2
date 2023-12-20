from __future__ import print_function

from core.CCP4PluginScript import CPluginScript
import sys, os
import pickle
from core import CCP4ErrorHandling
from core import CCP4Modules
from pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script import phaser_MR_AUTO
#from lxml import etree
from xml.etree import ElementTree as ET

class phaser_MR_FRF(phaser_MR_AUTO.phaser_MR_AUTO):

    TASKNAME = 'phaser_MR_FRF'                                  # Task name - should be same as class name
    TASKTITLE='Rotation function - PHASER'
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    ASYNCHRONOUS = False
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

        inputObject = phaser.InputMR_FRF()

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
        self.resultObject = phaser.runMR_FRF(inputObject, outputObject)
        self.finishCaptureCPlusPlusStdout()

        if not self.resultObject.Success():
            self.appendErrorReport(105, self.resultObject.ErrorName() + '-' + self.resultObject.ErrorMessage())
            return CPluginScript.FAILED
                
        if self.analyseResults(self.resultObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        resultObject = self.resultObject
        solutions = resultObject.getDotRlist()
        if len(solutions) > 0:
            picklePath = str(self.container.outputData.RFILEOUT.fullPath)
            if sys.version_info > (3,0):
                with open(picklePath,'wb') as pickleFile:
                    try:
                        pickle.dump(solutions, pickleFile)
                    except Exception as e:
                        print('Unable to Pickle Rfile solutions', e)
                    self.container.outputData.RFILEOUT.annotation.set('Rfile (pkl) from rotation search')
            else:
                with open(picklePath,'w') as pickleFile:
                    try:
                        pickle.dump(solutions, pickleFile)
                    except Exception as e:
                        print('Unable to Pickle Rfile solutions', e)
                    self.container.outputData.RFILEOUT.annotation.set('Rfile (pkl) from rotation search')
        #Remove old digested summaries and add new ones parsed from the result summary block
        parent_map = {c: p for p in self.xmlroot.iter() for c in p}
        for summaryNode in self.xmlroot.findall('Summary'):
            parent_map[summaryNode].remove(summaryNode)
        summary_buffer = '***'
        for text in resultObject.summary().split('\n'):
            if text.startswith("**********") and not summary_buffer.strip().endswith("***"):
                summaryNode = ET.SubElement(self.xmlroot,'Summary')
                summaryNode.text = summary_buffer
                summary_buffer = ""
            summary_buffer += (text+'\n')
        summaryNode = ET.SubElement(self.xmlroot,'Summary')
        summaryNode.text = summary_buffer
                
        self.flushXML(self.xmlroot)
        return CPluginScript.SUCCEEDED


    def analyseResults(self, results):
        solutionsNode = ET.SubElement(self.xmlroot,'PhaserMrSolutions')
        if len(results.getDotRlist()) > 0:
            node=self.subElementWithNameAndText(solutionsNode,'solutionsFound','True')
        else:
            node=self.subElementWithNameAndText(solutionsNode,'solutionsFound','False')
        solutionListNode = ET.SubElement(solutionsNode,'Solutions')

        for rlists in results.getDotRlist():
            for solution in rlists.RLIST:
                solutionNode = ET.SubElement(solutionListNode,'Solution')
                for property in ['DEEP', 'EULER', 'GRF', 'MODLID', 'RF', 'RFZ']:
                    value = getattr(solution,property,None)
                    if value is not None:
                        node = self.subElementWithNameAndText(solutionNode, property,str(value))

        self.flushXML(self.xmlroot)
        return  CPluginScript.SUCCEEDED

