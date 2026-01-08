from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.phaser_pipeline.script import phaser_pipeline


class phaser_rnp_pipeline(phaser_pipeline.phaser_pipeline):

    TASKNAME = 'phaser_rnp_pipeline'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    
    ERROR_CODES = {  200 : { 'description' : 'Phaser exited with error statut' }, 202 : { 'description' : 'Failed in harvest operation' }, 203 : { 'description' : 'Columns not present' }, 204 : { 'description' : 'Failed in plugin:',205 : { 'description' : 'Failed in pointless reindex operation' }, }, }
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)        
        self.checkOutputData()

        self.xmlroot = etree.Element('PhaserPipeline')

        self.F_SIGF_TOUSE = self.container.inputData.F_SIGF
        self.FREERFLAG_TOUSE = self.container.inputData.FREERFLAG
        rv = self.runPointless()
        #print 'self.F_SIGF_TOUSE is',self.F_SIGF_TOUSE
        rv = self.runPhaser(F_SIGF=self.F_SIGF_TOUSE)
        XYZIN_TOUSE = self.container.outputData.XYZOUT[0]

        if self.container.inputData.RUNREFMAC:
            self.runRefmac(F_SIGF=self.F_SIGF_TOUSE, FREERFLAG=self.FREERFLAG_TOUSE, XYZIN=XYZIN_TOUSE)

        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def createEnsembleElements(self):
        from ccp4i2.core.CCP4ModelData import CAtomSelection, CPdbDataFile, CPdbEnsembleItem
        elements = self.container.inputData.ENSEMBLES
        #Before removing all elements from this list, I have to set its listMinLength to 0
        self.container.inputData.ENSEMBLES.setQualifiers({'listMinLength':0})
        while len(elements) > 0: elements.remove(elements[-1])

        iEntity = 1
        if len(self.container.inputData.SELECTIONS) == 0:
            self.container.inputData.SELECTIONS.append(self.container.inputData.SELECTIONS.makeItem())
            self.container.inputData.SELECTIONS[-1].text.set("/*/*/*.*")
        for selection in self.container.inputData.SELECTIONS:
            self.container.inputData.ENSEMBLES.append(self.container.inputData.ENSEMBLES.makeItem())
            ensemble = self.container.inputData.ENSEMBLES[-1]
            ensemble.number.set(1)
            ensemble.label.set('Fragment_'+str(iEntity))
            elements = ensemble.pdbItemList
            while len(elements) > 1: elements.remove(elements[-1])
            pdbItem = elements[-1]
            pdbItem.structure.set(self.container.inputData.XYZIN_PARENT)
            pdbItem.structure.selection.text.set(str(selection))
            pdbItem.identity_to_target.set(0.9)
            self.container.inputData.USINGSOLELEMENTS.append(self.container.inputData.USINGSOLELEMENTS.makeItem())
            self.container.inputData.USINGSOLELEMENTS[-1].set('Fragment_'+str(iEntity))
            iEntity += 1

    def runPhaser(self, F_SIGF=None):
        #print 'F_SIGF is ',F_SIGF
        try:
            self.createEnsembleElements()
            phaserPlugin = self.makePluginObject('phaser_MR_RNP')
            #This funky arrangement is the way to ensure that the plugin behaves the same
            #when it is a part of the ipeline as it does when it is run alone...something about defaults I guess
            for attrName in phaserPlugin.container.keywords.dataOrder():
                if hasattr(self.container.keywords,attrName):
                    attr = getattr(self.container.keywords,attrName)
                    if hasattr(attr,'isSet') and attr.isSet():
                        setattr(phaserPlugin.container.keywords,attrName,attr)
            phaserPlugin.container.inputData=self.container.inputData
            phaserPlugin.container.inputData.RESOLUTION_LOW.set(25.0)
            phaserPlugin.container.inputData.RESOLUTION_HIGH.set(3.0)
            if F_SIGF is not None: phaserPlugin.container.inputData.F_SIGF=F_SIGF
            
            phaserPlugin.container.inputData.F_SIGF.loadFile()
            columns = phaserPlugin.container.inputData.F_SIGF.fileContent.getListOfColumns()
            columnStrings = [column.columnLabel.__str__() for column in columns]
            print(columnStrings)
            if 'F' in columnStrings: phaserPlugin.container.inputData.F_OR_I.set('F')

            rv = phaserPlugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(204,'phaser_MR_RNP')
                self.reportStatus(rv)
            
            pluginOutputs=phaserPlugin.container.outputData
            pipelineOutputs = self.container.outputData
            self.appendXML(phaserPlugin.makeFileName('PROGRAMXML'),'PhaserMrResults')
            self.harvestFile(pluginOutputs.SOLOUT, pipelineOutputs.SOLOUT)
            for outputListType in ['XYZOUT', 'MAPOUT', 'DIFMAPOUT','PHASEOUT']:
                pluginOutputList = getattr(pluginOutputs, outputListType, None)
                pipelineOutputList = getattr(pipelineOutputs, outputListType, None)
                self.harvestList(pluginOutputList, pipelineOutputList)
        except:
            self.appendErrorReport(202,'phaser_MR_RNP')
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED


    def runPointless(self):
        try:
            pointlessPlugin = self.makePluginObject('pointless_reindexToMatch')
            pointInp = pointlessPlugin.container.inputData
            pointlessPlugin.container.controlParameters.REFERENCE = 'XYZIN_REF'
            pointInp.XYZIN_REF = self.container.inputData.XYZIN_PARENT
            pointInp.F_SIGF.set(self.container.inputData.F_SIGF)
            if self.container.inputData.FREERFLAG.isSet():
                pointInp.FREERFLAG = self.container.inputData.FREERFLAG
            rv = pointlessPlugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(204,'pointless_reindexToMatch')
                self.reportStatus(rv)
            
            pluginOutputs = pointlessPlugin.container.outputData
            pipelineOutputs = self.container.outputData
            
            pluginOutputs.F_SIGF_OUT.loadFile()
            self.container.inputData.F_SIGF.loadFile()
            cellsAreSame = pluginOutputs.F_SIGF_OUT.fileContent.clipperSameCell(self.container.inputData.F_SIGF.fileContent)
            #print '\n\n\n ClipperSameCell',cellsAreSame,'\n'
            # pluginOutputs.F_SIGF_OUT.fileContent
            #print self.container.inputData.F_SIGF.fileContent
            if not cellsAreSame['validity']:
                #print 'self.F_SIGF_TOUSE was',self.F_SIGF_TOUSE
                self.harvestFile(pluginOutputs.F_SIGF_OUT, pipelineOutputs.F_SIGF_OUT)
                self.F_SIGF_TOUSE = pluginOutputs.F_SIGF_OUT
                #print 'self.F_SIGF_TOUSE became',self.F_SIGF_TOUSE
                #print 'pipelineOutputs.F_SIGF_OUT is',pipelineOutputs.F_SIGF_OUT
                if self.container.inputData.FREERFLAG.isSet():
                    self.harvestFile(pluginOutputs.FREERFLAG_OUT, pipelineOutputs.FREERFLAG_OUT)
                    self.FREERFLAG_TOUSE = pluginOutputs.FREERFLAG_OUT
            
            try:
                self.appendXML(pointlessPlugin.makeFileName('PROGRAMXML'),'Pointless')
            except:
                pass
        except:
            self.appendErrorReport(205,'pointless_reindexToMatch')
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED
