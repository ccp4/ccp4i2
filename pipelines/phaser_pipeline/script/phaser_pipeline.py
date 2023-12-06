from __future__ import print_function

try:
    import ccp4mg
    import mmdb2 as mmdb
except:
    print('FAILED CCP4ModelData imported ccp4mg')
import mmut

from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
import sys, os
from core import CCP4ErrorHandling
from core import CCP4Modules
from lxml import etree
from core import CCP4Utils
  
class phaser_pipeline(CPluginScript):

    TASKNAME = 'phaser_pipeline'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    ASYNCHRONOUS = False
    PERFORMANCECLASS = 'CRefinementPerformance'
    SEPARATEDATA=True
    INTERRUPTABLE=True

    ERROR_CODES = {  200 : { 'description' : 'Phaser exited with error status' }, 202 : { 'description' : 'Failed in harvest operation' }, 203 : { 'description' : 'Columns not present' }, 204 : { 'description' : 'Failed in plugin:' },}
    WHATNEXT = ['prosmart_refmac','buccaneer_build_refine_mr','coot_rebuild','modelcraft']
    

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    
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
        self.phaserPlugin.doAsync = False
        #self.phaserPlugin.waitForFinished = -1
        #self.phaserPlugin.setFinishHandler(self.phaser_MR_AUTO_Finished)
        self.connectSignal(self.phaserPlugin,'finished', self.phaserFinished)
        self.oldXMLLength = 0
        self.phaserPlugin.callbackObject.addResponder(self.phaserXMLUpdated)
        print('Off to see the wizard')
        rv = self.phaserPlugin.process()
        sys.stdout.flush()
        print('Rv', rv)
        sys.stdout.flush()
        if rv == CPluginScript.FAILED:
            print('Oh')
            with open(self.phaserPlugin.makeFileName('LOG'),"r") as logFile:
                wasInterrupted = 'KILL-FILE DETECTED ERROR' in logFile.read()
                print('Was interrupted',wasInterrupted)
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
        self.renameFile(tmpFilename,self.makeFileName('PROGRAMXML'))

    @QtCore.Slot(dict)
    def phaserFinished(self, statusDict = {}):
        sys.stdout.flush()
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
        print('StatusDict',statusDict)
        sys.stdout.flush()
        self.checkSolutionsFound(statusDict=statusDict, failedErrCode=200)
        if len(self.phaserPlugin.container.outputData.XYZOUT) > 0:
            self.checkFinishStatus(statusDict=statusDict,failedErrCode=200,outputFile = self.phaserPlugin.container.outputData.XYZOUT[0] ,noFileErrCode=207)
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

        if self.container.inputData.RUNCOOT:
            self.runCoot(MAPIN=self.container.outputData.MAPOUT[0], XYZIN=XYZIN_TOUSE)
            XYZIN_TOUSE = self.container.outputData.XYZOUT_COOT
        
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
        except:
            self.appendErrorReport(204,'pointless_reindexToMatch')
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED

    def runCsymmatch(self):
        try:
            csymmatchPlugin = self.makePluginObject('csymmatch')
            csymmatchInp.set(csymmatchPlugin.container.inputData)
            csymmatchInp.XYZIN_QUERY.set(self.container.outputData.XYZOUT[0])
            csymmatchInp.XYZIN_TARGET.set(self.container.inputData.XYZIN_TARGET)
            rv = csymmatchPlugin.process()
            if rv != CPluginScript.SUCCEEDED: self.reportStatus(rv)
            
            pluginOutputs = csymmatchPlugin.container.outputData
            pipelineOutputs = self.container.outputData

            self.harvestFile(pluginOutputs.XYZOUT, pipelineOutputs.XYZOUT_CSYMMATCH)
            self.appendXML(csymmatchPlugin.makeFileName('PROGRAMXML'),'Csymmatch')
            pipelineOutputs.XYZOUT_CSYMMATCH.annotation.set('Coordinates moved to match reference structure')
        except:
            self.appendErrorReport(204,'csymmatch')
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED

    def runCoot(self, MAPIN=None, XYZIN=None):
        try:
            try:
                cootPlugin = self.makePluginObject('coot_script_lines')
                xyzinList = cootPlugin.container.inputData.XYZIN
                xyzinList.append(xyzinList.makeItem())
                xyzinList[-1].set(XYZIN)
                fphiinList = cootPlugin.container.inputData.FPHIIN
                fphiinList.append(fphiinList.makeItem())
                fphiinList[-1].set(MAPIN)
                cootPlugin.container.controlParameters.SCRIPT = '''fill_partial_residues(MolHandle_1)
fit_protein(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
            except:
                self.appendErrorReport(204,'coot_script_lines -setup')
                self.reportStatus(CPluginScript.FAILED)
            try:
                cootPlugin.doAsync=False
                rv = cootPlugin.process()
                if rv != CPluginScript.SUCCEEDED: self.reportStatus(rv)
            except:
                self.appendErrorReport(204,'coot_script_lines -execute')
                self.reportStatus(CPluginScript.FAILED)
            try:
                pluginOutputs = cootPlugin.container.outputData
                pipelineOutputs = self.container.outputData

                self.harvestFile(pluginOutputs.XYZOUT[0], pipelineOutputs.XYZOUT_COOT)
                self.appendXML(cootPlugin.makeFileName('PROGRAMXML'),'coot_script_lines')
                pipelineOutputs.XYZOUT_COOT.annotation.set('Coordinates filled and fitted by COOT')
            except:
                self.appendErrorReport(204,'coot_script_lines -postprocess')
                self.reportStatus(CPluginScript.FAILED)
        except:
            self.appendErrorReport(204,'coot_script_lines')
            self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.SUCCEEDED

    def runSheetbend(self, F_SIGF=None, FREERFLAG=None, XYZIN=None):
        try:
            self.sheetbendPlugin = self.makePluginObject('sheetbend')
            if XYZIN is not None: self.sheetbendPlugin.container.inputData.XYZIN.set(XYZIN)
            if F_SIGF is not None: self.sheetbendPlugin.container.inputData.F_SIGF.set(F_SIGF)
            if FREERFLAG is not None: self.sheetbendPlugin.container.inputData.FREERFLAG.set(FREERFLAG)
            self.sheetbendPlugin.doAsync=False
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
            self.refmacPlugin.doAsync = False
            rv = self.refmacPlugin.process()
            if rv == CPluginScript.FAILED: self.reportStatus(rv)
            
            pluginOutputs=self.refmacPlugin.container.outputData
            pipelineOutputs = self.container.outputData
            self.harvestFile(pluginOutputs.FPHIOUT, pipelineOutputs.MAPOUT_REFMAC)
            self.harvestFile(pluginOutputs.DIFFPHIOUT, pipelineOutputs.DIFMAPOUT_REFMAC)
            self.harvestFile(pluginOutputs.XYZOUT, pipelineOutputs.XYZOUT_REFMAC)
            
            self.appendXML(self.refmacPlugin.makeFileName('PROGRAMXML'),'REFMAC')
        except:
            self.appendErrorReport(204,'refmac Preamble')
            self.reportStatus(CPluginScript.FAILED)
        try:
            self.container.outputData.PERFORMANCEINDICATOR.set(self.refmacPlugin.container.outputData.PERFORMANCEINDICATOR)
        except:
            self.appendErrorReport(204,'refmac PI')
            self.reportStatus(CPluginScript.FAILED)
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
        except:
            self.appendErrorReport(202,str(pluginOutputItem.fullPath)+' '+str(pipelineOutputItem.fullPath))
            self.reportStatus(CPluginScript.FAILED)

    def appendXML(self, changedFile, replacingElementOfType=None):
        for oldNode in self.xmlroot.xpath(replacingElementOfType):
            self.xmlroot.remove(oldNode)
        try:
            newXML = CCP4Utils.openFileToEtree(changedFile)
        except:
            newXML = etree.Element(replacingElementOfType)
        self.xmlroot.append(newXML)
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlfile:
            CCP4Utils.writeXML(xmlfile,etree.tostring(self.xmlroot,pretty_print=True))

    def checkFinishStatus( self, statusDict,failedErrCode,outputFile = None,noFileErrCode= None):
        if len(statusDict)>0 and statusDict['finishStatus'] == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(statusDict['finishStatus'])
        try:
            assert outputFile.exists(),'Entity provided is not CDataFile or does not exist'
        except:
            self.appendErrorReport(noFileErrCode,'Expected file: '+str(outputFile))
            self.reportStatus(CPluginScript.FAILED)

    def checkSolutionsFound(self, statusDict, failedErrCode):
        if len(statusDict)>0 and statusDict['finishStatus'] == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(statusDict['finishStatus'])
        self.appendXML(self.phaserPlugin.makeFileName('PROGRAMXML'),'PhaserMrResults')
        if self.xmlroot.xpath('//solutionsFound')[0].text == 'False':
            self.reportStatus(CPluginScript.UNSATISFACTORY)
