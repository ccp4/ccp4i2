from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script import phaser_MR_AUTO


class phaser_MR_RNP(phaser_MR_AUTO.phaser_MR_AUTO):

    TASKNAME = 'phaser_MR_RNP'                                  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' },}

    def startProcess(self):
        import phaser
        outputObject = phaser.Output()
        outputObject.setPhenixCallback(self.callbackObject)

        resultObject = self.runMR_DAT(outputObject)
        if resultObject == CPluginScript.FAILED: return CPluginScript.FAILED
        
        self.inputHall = resultObject.getSpaceGroupHall()
        self.inputSpaceGroup = resultObject.getSpaceGroupName()

        inputObject = phaser.InputMR_RNP()
        
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

        if self.parseSolutions(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        
        self.prepareCaptureCPlusPlusStdoutToLog()
        inputObject.setMUTE(False)
        self.resultObject = phaser.runMR_RNP(inputObject, outputObject)
        self.finishCaptureCPlusPlusStdout()

        if not self.resultObject.Success():
            self.appendErrorReport(105, self.resultObject.ErrorName() + '-' + self.resultObject.ErrorMessage())
            return CPluginScript.FAILED
        
        if self.analyseResults(self.resultObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        
        return CPluginScript.SUCCEEDED

    def parseSolutions(self, inputObject):
        if self.container.inputData.SOLIN.isSet():
            return super(phaser_MR_RNP,self).parseSolutions(inputObject)
        else:
            inputObject.addSOLU_SET('Fragments')
            for usingSolution in self.container.inputData.USINGSOLELEMENTS:
                fragName = str(usingSolution)
                print('ghb', fragName)
                inputObject.addSOLU_6DIM_ENSE(fragName, [0.,0.,0.],False,[0.,0.,0.],20.0, False, False, False, 1,1)
            return CPluginScript.SUCCEEDED
