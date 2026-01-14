import os

from lxml import etree

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class phaser_pipeline(CPluginScript):

    TASKNAME = 'phaser_pipeline'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    PERFORMANCECLASS = 'CRefinementPerformance'
    SEPARATEDATA=True
    INTERRUPTABLE=True

    ERROR_CODES = {
        200: {'description': 'Phaser exited with error status'},
        201: {'description': 'Exception in phaser pipeline setup'},
        202: {'description': 'Failed in harvest operation'},
        203: {'description': 'Columns not present'},
        204: {'description': 'Failed in plugin'},
        205: {'description': 'Exception in pointless_reindexToMatch'},
        206: {'description': 'Exception in csymmatch'},
        207: {'description': 'No output files in list'},
        209: {'description': 'Exception in sheetbend'},
        210: {'description': 'Exception in refmac'},
        211: {'description': 'Exception in harvestFile'},
    }
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()
        self.xmlroot = etree.Element('PhaserPipeline')
        if self.container.inputData.MODE_TY == "MR_FRF":
            self.phaserPlugin = self.makePluginObject('phaser_MR_FRF')
        elif self.container.inputData.MODE_TY == "MR_FTF":
            self.phaserPlugin = self.makePluginObject('phaser_MR_FTF')
        elif self.container.inputData.MODE_TY == "MR_PAK":
            self.phaserPlugin = self.makePluginObject('phaser_MR_PAK')
        elif self.container.inputData.MODE_TY == "MR_RNP":
            self.phaserPlugin = self.makePluginObject('phaser_MR_RNP')
        else:
            self.phaserPlugin = self.makePluginObject('phaser_MR_AUTO')

        #This funky arrangement is the way to ensure that the plugin behaves the same
        #when it is a part of the pipeline as it does when it is run alone...something about defaults I guess
        #Note without allSet=False isSet() returns False for a CContainer containing any items that are unset
        for attrName in self.phaserPlugin.container.keywords.dataOrder():
            if hasattr(self.container.keywords,attrName):
                attr = getattr(self.container.keywords,attrName)
                if hasattr(attr,'isSet') and attr.isSet(allSet=False):
                    print("Setting",attrName,attr)
                    #setattr(self.phaserPlugin.container.keywords,attrName,attr)
                    getattr(self.phaserPlugin.container.keywords,attrName).set(attr)
        self.phaserPlugin.container.inputData.set(self.container.inputData)
        self.phaserPlugin.container.inputData.KILLFILEPATH.set(os.path.join(self.getWorkDirectory(),'INTERRUPT'))
        self.oldXMLLength = 0
        self.phaserPlugin.callbackObject.addResponder(self.phaserXMLUpdated)
        rv = self.phaserPlugin.process()
        self.phaserFinished(rv)
        if rv == CPluginScript.FAILED:
            # Check if LOG file exists before reading it
            # In standalone/i2run mode, subjobs with RUNEXTERNALPROCESS=False don't create LOG files
            # In GUI mode, subjobs run as proper database jobs and do create LOG files
            log_file_path = self.phaserPlugin.makeFileName('LOG')
            wasInterrupted = False
            if os.path.exists(log_file_path):
                with open(log_file_path, "r") as logFile:
                    wasInterrupted = 'KILL-FILE DETECTED ERROR' in logFile.read()
            if wasInterrupted:
                self.reportStatus('CPluginScript.INTERRUPTED')
                return CPluginScript.INTERRUPTED
            else:
                self.reportStatus(rv)
                return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def phaserXMLUpdated(self, newXML):
        for oldNode in self.xmlroot.xpath('PhaserMrResults'): self.xmlroot.remove(oldNode)
        from copy import deepcopy
        self.xmlroot.append(deepcopy(newXML))
        tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
        with open(tmpFilename,'w') as xmlfile:
            CCP4Utils.writeXML(xmlfile,etree.tostring(self.xmlroot,pretty_print=True))
        finalFilename = self.makeFileName('PROGRAMXML')
        self.renameFile(tmpFilename,finalFilename)

    def phaserFinished(self, finishStatus):
        # If phaser subjob failed, propagate the failure status
        if finishStatus == CPluginScript.FAILED:
            self.reportStatus(CPluginScript.FAILED)
            return

        if self.container.inputData.MODE_TY in ['MR_FRF', 'MR_FTF', 'MR_PAK']:
            pluginOutputs=self.phaserPlugin.container.outputData
            pipelineOutputs = self.container.outputData
            if self.container.inputData.MODE_TY == 'MR_FRF':
                self.harvestFile(pluginOutputs.RFILEOUT, pipelineOutputs.RFILEOUT)
            else:
                self.harvestFile(pluginOutputs.SOLOUT, pipelineOutputs.SOLOUT)
            self.appendXML(self.phaserPlugin.makeFileName('PROGRAMXML'),'PhaserMrResults')
            self.reportStatus(CPluginScript.SUCCEEDED)
            return
        self.checkSolutionsFound(finishStatus=finishStatus, failedErrCode=200)
        if len(self.phaserPlugin.container.outputData.XYZOUT) > 0:
            self.checkFinishStatus(finishStatus=finishStatus,failedErrCode=200,outputFile = self.phaserPlugin.container.outputData.XYZOUT[0] ,noFileErrCode=207)
        else:
            self.appendErrorReport(207,'No output files in list')
            self.reportStatus(CPluginScript.FAILED)
        self.harvestPhaserPlugin()
        self.appendXML(self.phaserPlugin.makeFileName('PROGRAMXML'),'PhaserMrResults')
        F_SIGF_TOUSE = self.container.inputData.F_SIGF
        FREERFLAG_TOUSE = self.container.inputData.FREERFLAG
        XYZIN_TOUSE = self.container.outputData.XYZOUT[0]

        if self.phaserPlugin.container.outputData.dataReindexed:
            self.runPointless()
            F_SIGF_TOUSE = self.container.outputData.F_SIGF_OUT
            FREERFLAG_TOUSE = self.container.outputData.FREERFLAG_OUT

        if self.container.inputData.XYZIN_TARGET.isSet():
            self.runCsymmatch()
            XYZIN_TOUSE = self.container.outputData.XYZOUT_CSYMMATCH

        if self.container.inputData.RUNSHEETBEND:
            self.runSheetbend(F_SIGF=F_SIGF_TOUSE, FREERFLAG=FREERFLAG_TOUSE, XYZIN=XYZIN_TOUSE)
            XYZIN_TOUSE = self.container.outputData.XYZOUT_SHEETBEND

        if self.container.inputData.RUNREFMAC:
            self.runRefmac(F_SIGF=F_SIGF_TOUSE, FREERFLAG=FREERFLAG_TOUSE, XYZIN=XYZIN_TOUSE)

        self.reportStatus(CPluginScript.SUCCEEDED)

    def harvestPhaserPlugin(self):
        pluginOutputs=self.phaserPlugin.container.outputData
        pipelineOutputs = self.container.outputData
        self.harvestFile(pluginOutputs.SOLOUT, pipelineOutputs.SOLOUT)
        for outputListType in ['XYZOUT', 'MAPOUT', 'DIFMAPOUT','PHASEOUT']:
            pluginOutputList = getattr(pluginOutputs, outputListType, None)
            pipelineOutputList = getattr(pipelineOutputs, outputListType, None)
            self.harvestList(pluginOutputList, pipelineOutputList)

    def runPointless(self):
        try:
            pointlessPlugin = self.makePluginObject('pointless_reindexToMatch')
            pointInp = pointlessPlugin.container.inputData
            pointlessPlugin.container.controlParameters.REFERENCE = 'HKLIN_FMAP_REF'
            pointInp.HKLIN_FMAP_REF.set(self.container.outputData.MAPOUT[0])
            pointInp.F_SIGF.set(self.container.inputData.F_SIGF)
            if self.container.inputData.FREERFLAG.isSet():
                pointInp.FREERFLAG.set(self.container.inputData.FREERFLAG)
            rv = pointlessPlugin.process()
            if rv != CPluginScript.SUCCEEDED: self.reportStatus(rv)

            pluginOutputs = pointlessPlugin.container.outputData
            pipelineOutputs = self.container.outputData

            self.harvestFile(pluginOutputs.F_SIGF_OUT, pipelineOutputs.F_SIGF_OUT)
            if self.container.inputData.FREERFLAG.isSet():
                self.harvestFile(pluginOutputs.FREERFLAG_OUT, pipelineOutputs.FREERFLAG_OUT)
        except Exception as e:
            self.appendErrorReport(205, 'Exception in pointless_reindexToMatch: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED

    def runCsymmatch(self):
        try:
            csymmatchPlugin = self.makePluginObject('csymmatch')
            csymmatchInp = csymmatchPlugin.container.inputData
            csymmatchInp.XYZIN_QUERY.set(self.container.outputData.XYZOUT[0])
            csymmatchInp.XYZIN_TARGET.set(self.container.inputData.XYZIN_TARGET)
            rv = csymmatchPlugin.process()
            if rv != CPluginScript.SUCCEEDED: self.reportStatus(rv)

            pluginOutputs = csymmatchPlugin.container.outputData
            pipelineOutputs = self.container.outputData

            self.harvestFile(pluginOutputs.XYZOUT, pipelineOutputs.XYZOUT_CSYMMATCH)
            self.appendXML(csymmatchPlugin.makeFileName('PROGRAMXML'),'Csymmatch')
            pipelineOutputs.XYZOUT_CSYMMATCH.annotation.set('Coordinates moved to match reference structure')
        except Exception as e:
            self.appendErrorReport(206, 'Exception in csymmatch: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED

    def runSheetbend(self, F_SIGF=None, FREERFLAG=None, XYZIN=None):
        try:
            self.sheetbendPlugin = self.makePluginObject('sheetbend')
            if XYZIN is not None: self.sheetbendPlugin.container.inputData.XYZIN.set(XYZIN)
            if F_SIGF is not None: self.sheetbendPlugin.container.inputData.F_SIGF.set(F_SIGF)
            if FREERFLAG is not None: self.sheetbendPlugin.container.inputData.FREERFLAG.set(FREERFLAG)
            rv = self.sheetbendPlugin.process()
            if rv == CPluginScript.FAILED: self.reportStatus(rv)
            pluginOutputs=self.sheetbendPlugin.container.outputData
            pipelineOutputs = self.container.outputData
            self.harvestFile(pluginOutputs.XYZOUT, pipelineOutputs.XYZOUT_SHEETBEND)
            pipelineOutputs.XYZOUT_SHEETBEND.annotation.set('Atomic model after Shift field refinement')
            self.appendXML(self.sheetbendPlugin.makeFileName('PROGRAMXML'), 'SheetbendResult')
        except Exception as e:
            self.appendErrorReport(204,'sheetbend: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED

    def runRefmac(self,F_SIGF=None,FREERFLAG=None,XYZIN=None):
        try:
            # refmac wrapper run with 10 cycles
            self.refmacPlugin = self.makePluginObject('refmac')
            if XYZIN is not None: self.refmacPlugin.container.inputData.XYZIN.set(XYZIN)
            if F_SIGF is not None: self.refmacPlugin.container.inputData.F_SIGF.set(F_SIGF)
            if FREERFLAG is not None: self.refmacPlugin.container.inputData.FREERFLAG.set(FREERFLAG)
            self.refmacPlugin.container.controlParameters.HYDROGENS.set('NO')
            self.refmacPlugin.container.controlParameters.NCYCLES.set(10)
            self.refmacPlugin.container.controlParameters.PHOUT.set(False)
            self.refmacPlugin.container.controlParameters.USE_JELLY.set(True)
            self.refmacPlugin.container.controlParameters.JELLY_SIGMA.set(0.05)
            self.refmacPlugin.container.controlParameters.MAKE_NEW_LIGAND_EXIT.set(False)
            rv = self.refmacPlugin.process()
            if rv == CPluginScript.FAILED: self.reportStatus(rv)

            pluginOutputs=self.refmacPlugin.container.outputData
            pipelineOutputs = self.container.outputData
            self.harvestFile(pluginOutputs.FPHIOUT, pipelineOutputs.MAPOUT_REFMAC)
            self.harvestFile(pluginOutputs.DIFFPHIOUT, pipelineOutputs.DIFMAPOUT_REFMAC)
            self.harvestFile(pluginOutputs.XYZOUT, pipelineOutputs.XYZOUT_REFMAC)

            self.appendXML(self.refmacPlugin.makeFileName('PROGRAMXML'),'REFMAC')
        except Exception as e:
            self.appendErrorReport(210, 'Exception in refmac: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        try:
            self.container.outputData.PERFORMANCEINDICATOR.set(self.refmacPlugin.container.outputData.PERFORMANCEINDICATOR)
        except Exception as e:
            self.appendErrorReport(210, 'Exception copying refmac performance indicator: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def harvestList(self, pluginOutputList, pipelineOutputList):
        for pluginOutputListItem in pluginOutputList:
            pipelineOutputList.append(pipelineOutputList.makeItem())
            pipelineOutputListItem = pipelineOutputList[-1]
            pipelineOutputListItem.fullPath = os.path.join(self.workDirectory,os.path.basename(str(pluginOutputListItem.fullPath)))
            self.harvestFile(pluginOutputListItem, pipelineOutputListItem)

    def harvestFile(self, pluginOutputItem, pipelineOutputItem):
        import shutil
        try:
            shutil.copyfile(str(pluginOutputItem.fullPath), str(pipelineOutputItem.fullPath))
            pipelineOutputItem.annotation.set(pluginOutputItem.annotation)
            pipelineOutputItem.contentFlag.set(pluginOutputItem.contentFlag)
            pipelineOutputItem.subType.set(pluginOutputItem.subType)
        except Exception as e:
            self.appendErrorReport(211, 'Exception in harvestFile: ' + str(pluginOutputItem.fullPath) + ' -> ' + str(pipelineOutputItem.fullPath) + ': ' + str(e))
            self.reportStatus(CPluginScript.FAILED)

    def appendXML(self, changedFile, replacingElementOfType=None):
        for oldNode in self.xmlroot.xpath(replacingElementOfType):
            self.xmlroot.remove(oldNode)
        try:
            newXML = CCP4Utils.openFileToEtree(changedFile)
        except Exception as e:
            newXML = etree.Element(replacingElementOfType)
        self.xmlroot.append(newXML)
        output_file = self.makeFileName('PROGRAMXML')
        with open(output_file,'w') as xmlfile:
            CCP4Utils.writeXML(xmlfile,etree.tostring(self.xmlroot,pretty_print=True))

    def checkFinishStatus( self, finishStatus,failedErrCode,outputFile = None,noFileErrCode= None):
        if finishStatus == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(finishStatus)
        try:
            assert outputFile.exists(),'Entity provided is not CDataFile or does not exist'
        except Exception as e:
            self.appendErrorReport(noFileErrCode,'Expected file: '+str(outputFile) + ' - ' + str(e))
            self.reportStatus(CPluginScript.FAILED)

    def checkSolutionsFound(self, finishStatus, failedErrCode):
        if finishStatus == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(finishStatus)
        self.appendXML(self.phaserPlugin.makeFileName('PROGRAMXML'),'PhaserMrResults')
        if self.xmlroot.xpath('//solutionsFound')[0].text == 'False':
            self.reportStatus(CPluginScript.UNSATISFACTORY)
