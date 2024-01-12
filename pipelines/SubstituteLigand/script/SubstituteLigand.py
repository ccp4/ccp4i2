from __future__ import print_function

from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
from lxml import etree
import os
import shutil

class SubstituteLigand(CPluginScript):
    TASKNAME = 'SubstituteLigand'            # Task name - should be same as class name
    TASKVERSION= 0.0                    # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    WHATNEXT = ['coot_rebuild']
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    ERROR_CODES = { 201 : { 'description' : 'Failed in SubstituteLigand' },}

    def __init__(self, *args,**kws):
        super(SubstituteLigand, self).__init__(*args, **kws)
        self.xmlroot = ET.Element('SubstituteLigand')
        self.obsToUse = None
        self.freerToUse = None
    
        if self.container.controlParameters.OBSAS.__str__() != 'UNMERGED':
            #remove any (potentially invalid) entries from UNMERGED list
            while len(self.container.inputData.UNMERGEDFILES)>0:
                self.container.inputData.UNMERGEDFILES.remove(self.container.inputData.UNMERGEDFILES[-1])

    def process(self):
        invalidFiles = self.checkInputData()
        for invalidFile in invalidFiles:
            if self.container.controlParameters.OBSAS.__str__() == 'MERGED' and invalidFile.__str__() == 'UNMERGEDFILES':
                invalidFiles.remove(invalidFile)
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()

        # Chop out the chunk of file we want to use
        selAtomsFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected_atoms.pdb'))
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(selAtomsFilePath)
        from core.CCP4ModelData import CPdbDataFile
        self.selAtomsFile = CPdbDataFile(selAtomsFilePath)
        
        if self.container.controlParameters.LIGANDAS.__str__() == 'DICT':
            self.dictToUse = self.container.inputData.DICTIN
            self.dictDone()
        elif self.container.controlParameters.LIGANDAS.__str__() == 'NONE':
            self.dictDone()
        elif self.container.controlParameters.LIGANDAS.__str__() != 'NONE':
            self.lidiaAcedrgPlugin = self.makePluginObject('LidiaAcedrgNew')
            self.lidiaAcedrgPlugin.container.inputData.MOLSMILESORSKETCH = self.container.controlParameters.LIGANDAS
            if self.container.inputData.MOLIN.isSet():
                self.lidiaAcedrgPlugin.container.inputData.MOLIN=self.container.inputData.MOLIN
            if self.container.inputData.SMILESIN.isSet():
                self.lidiaAcedrgPlugin.container.inputData.SMILESIN=self.container.inputData.SMILESIN
            self.lidiaAcedrgPlugin.container.inputData.CONFORMERSFROM = 'RDKIT'
            self.lidiaAcedrgPlugin.container.inputData.TLC=self.container.inputData.TLC
            self.connectSignal(self.lidiaAcedrgPlugin,'finished',self.lidiaAcedrg_finished)
            self.lidiaAcedrgPlugin.process()

    @QtCore.Slot(dict)
    def lidiaAcedrg_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(CPluginScript.FAILED)
        pluginRoot = CCP4Utils.openFileToEtree(self.lidiaAcedrgPlugin.makeFileName('PROGRAMXML'))
        self.xmlroot.append(pluginRoot)
        self.flushXML()
        self.harvestFile(self.lidiaAcedrgPlugin.container.outputData.DICTOUT_LIST[0], self.container.outputData.DICTOUT)
        self.dictToUse = self.container.outputData.DICTOUT
        self.dictDone()
            
    def dictDone(self):
        if self.container.controlParameters.OBSAS.__str__() == 'UNMERGED':
            self.aimlessCycle = 0
            self.aimlessPipe()
        else:
            self.obsToUse = self.container.inputData.F_SIGF_IN
            self.freerToUse = self.container.inputData.FREERFLAG_IN
            self.rigidBodyPipeline()

    def aimlessPipe(self):
        self.aimlessPlugin = self.makePluginObject('aimless_pipe')
        self.aimlessPlugin.container.controlParameters.MODE.set('MATCH')
        self.aimlessPlugin.container.controlParameters.RESOLUTION_RANGE = self.container.controlParameters.RESOLUTION_RANGE
        self.aimlessPlugin.container.controlParameters.SCALING_PROTOCOL.set('DEFAULT')
        self.aimlessPlugin.container.controlParameters.ONLYMERGE.set(False)
        self.aimlessPlugin.container.controlParameters.REFERENCE_DATASET.set('XYZ')
        self.aimlessPlugin.container.inputData.copyData(self.container.inputData,['UNMERGEDFILES'])
        self.aimlessPlugin.container.inputData.XYZIN_REF = self.container.inputData.XYZIN
        self.aimlessPlugin.container.controlParameters.TOLERANCE.set(10.)
        if self.container.inputData.FREERFLAG_IN.isSet():
            self.aimlessPlugin.container.inputData.FREERFLAG = self.container.inputData.FREERFLAG_IN
        self.connectSignal(self.aimlessPlugin,'finished',self.aimlessPlugin_finished)
        self.aimlessPlugin.process()
        
    @QtCore.Slot(dict)
    def aimlessPlugin_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(CPluginScript.FAILED)
        
        pluginRoot = CCP4Utils.openFileToEtree(self.aimlessPlugin.makeFileName('PROGRAMXML'))
        self.xmlroot.append(pluginRoot)
        #Here check on the resolution estimate, and cut dataset back accordingly
        if self.aimlessCycle == 0:
            datasetresultnodes = pluginRoot.findall(".//Result/Dataset")
            if len(datasetresultnodes) == 0:
                self.appendErrorReport(201,'No result nodes found')
                self.reportStatus(CPluginScript.FAILED)
            datasetresultnode = datasetresultnodes[0]
            dataresonodes = datasetresultnode.findall("ResolutionHigh/Overall")
            reslimitnodes = datasetresultnode.findall("ResolutionLimitEstimate")
            if len(dataresonodes)<1 or len(reslimitnodes)< 1:
                self.appendErrorReport(201,'Unable to identify resolution estimate limits')
                self.reportStatus(CPluginScript.FAILED)
            datareso = datasetresultnode.findall("ResolutionHigh/Overall")[0].text
            reslimitnodes = datasetresultnode.findall("ResolutionLimitEstimate")
            if reslimitnodes[0].get("type") == "CChalf":
                reslimit = reslimitnodes[0].findall("MaximumResolution")[0].text
            else:
                self.appendErrorReport(201,'Reso detection failed (first estimate not on CCHalf)')
                self.reportStatus(CPluginScript.FAILED)
            if float(reslimit) == float(datareso):
                self.aimlessCyclesFinished()
            else:
                self.aimlessCycle = 1
                self.container.controlParameters.RESOLUTION_RANGE.end = float(reslimit)
                self.aimlessPipe()
        else:
            self.aimlessCyclesFinished()

    def aimlessCyclesFinished(self):
        self.flushXML()
        self.harvestFile(self.aimlessPlugin.container.outputData.FREEROUT, self.container.outputData.FREERFLAG_OUT)
        self.harvestFile(self.aimlessPlugin.container.outputData.HKLOUT[0], self.container.outputData.F_SIGF_OUT)
        self.obsToUse = self.container.outputData.F_SIGF_OUT
        self.freerToUse = self.container.outputData.FREERFLAG_OUT
        self.rigidBodyPipeline()

    def rigidBodyPipeline(self):
        inp = self.container.inputData
        if hasattr(inp,'PIPELINE') and inp.PIPELINE.isSet() and inp.PIPELINE.__str__() == 'DIMPLE':
            self.i2Dimple()
        else:
            self.phaser_rnp_pipeline()

    def phaser_rnp_pipeline(self):
        self.rnpPlugin = self.makePluginObject('phaser_rnp_pipeline')
        self.rnpPlugin.container.inputData.XYZIN_PARENT = self.selAtomsFile
        self.rnpPlugin.container.inputData.F_SIGF = self.obsToUse
        self.rnpPlugin.container.inputData.FREERFLAG = self.freerToUse
        self.rnpPlugin.container.inputData.SELECTIONS.append({'text':'/*/*/*/*','pdbFileKey':'XYZIN_PARENT'})
        self.connectSignal(self.rnpPlugin,'finished',self.rnpPlugin_finished)
        self.rnpPlugin.process()

    @QtCore.Slot(dict)
    def rnpPlugin_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(status)
        try:
            pluginRoot = CCP4Utils.openFileToEtree(self.rnpPlugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
            self.flushXML()
            self.harvestFile(self.rnpPlugin.container.outputData.MAPOUT_REFMAC, self.container.outputData.FPHIOUT)
            self.mapToUse = self.container.outputData.FPHIOUT
            self.harvestFile(self.rnpPlugin.container.outputData.DIFMAPOUT_REFMAC, self.container.outputData.DIFFPHIOUT)
            if os.path.isfile(str(self.rnpPlugin.container.outputData.FREERFLAG_OUT)):
                self.harvestFile(self.rnpPlugin.container.outputData.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)
                self.freerToUse = self.container.outputData.FREERFLAG_OUT
            if os.path.isfile(str(self.rnpPlugin.container.outputData.F_SIGF_OUT)):
                self.harvestFile(self.rnpPlugin.container.outputData.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
                self.obsToUse = self.container.outputData.F_SIGF_OUT
            if self.container.controlParameters.LIGANDAS.__str__() == 'NONE':
                self.harvestFile(self.rnpPlugin.container.outputData.XYZOUT[0], self.container.outputData.XYZOUT)
                self.reportStatus(CPluginScript.SUCCEEDED)
            else:
                self.coordinatesForCoot = self.rnpPlugin.container.outputData.XYZOUT_REFMAC
                self.cootAddLigand()
        except:
            self.appendErrorReport(201,'Failed in rnpplugin finished')
            self.reportStatus(CPluginScript.FAILED)
        
    def i2Dimple(self):
        self.i2Dimple = self.makePluginObject('i2Dimple')
        self.i2Dimple.container.inputData.XYZIN = self.selAtomsFile
        self.i2Dimple.container.inputData.F_SIGF = self.obsToUse
        self.i2Dimple.container.inputData.FREERFLAG = self.freerToUse
        self.connectSignal(self.i2Dimple,'finished',self.i2Dimple_finished)
        self.i2Dimple.process()

    @QtCore.Slot(dict)
    def i2Dimple_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(status)
        try:
            pluginRoot = CCP4Utils.openFileToEtree(self.i2Dimple.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
            self.flushXML()
            self.harvestFile(self.i2Dimple.container.outputData.FPHIOUT, self.container.outputData.FPHIOUT)
            self.mapToUse = self.container.outputData.FPHIOUT
            self.harvestFile(self.i2Dimple.container.outputData.DIFFPHIOUT, self.container.outputData.DIFFPHIOUT)
            #Adopt the reindexed output of the dimple file if present
            if os.path.isfile(self.i2Dimple.container.outputData.F_SIGF_OUT.fullPath.__str__()):
                self.harvestFile(self.i2Dimple.container.outputData.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
            if os.path.isfile(self.i2Dimple.container.outputData.FREERFLAG_OUT.fullPath.__str__()):
                self.harvestFile(self.i2Dimple.container.outputData.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)
            if self.container.controlParameters.LIGANDAS.__str__() == 'NONE':
                self.harvestFile(self.i2Dimple.container.outputData.XYZOUT, self.container.outputData.XYZOUT)
                self.reportStatus(CPluginScript.SUCCEEDED)
            else:
                self.coordinatesForCoot = self.i2Dimple.container.outputData.XYZOUT
                self.cootAddLigand()
        except:
            self.appendErrorReport(201,'Failed in i2Dimple finished')
            self.reportStatus(CPluginScript.FAILED)
        
    def cootAddLigand(self):
        self.cootPlugin = self.makePluginObject('coot_script_lines')
        xyzinList = self.cootPlugin.container.inputData.XYZIN
        xyzinList.append(xyzinList.makeItem())
        xyzinList[-1].set(self.coordinatesForCoot)
        fphiinList = self.cootPlugin.container.inputData.FPHIIN
        fphiinList.append(fphiinList.makeItem())
        fphiinList[-1].set(self.mapToUse)
        self.cootPlugin.container.inputData.DICT = self.dictToUse
        #coot_stepped_refine,coot_fit_residues,coot_script_lines
        self.cootPlugin.container.controlParameters.SCRIPT = '''#Script to fit lignad into density
monomerMolNo = coot.get_monomer('LIG')
coot.add_ligand_clear_ligands()
coot.set_ligand_search_protein_molecule(MolHandle_1)
coot.set_ligand_search_map_molecule(MapHandle_1)
coot.add_ligand_search_wiggly_ligand_molecule(monomerMolNo)
#Execute search
nToCopy = 0
ligandsFound=coot.execute_ligand_search()
if ligandsFound is not False:
    nToCopy = len(ligandsFound)
#Check on ncs
equivs = coot.ncs_chain_ids(0)
if equivs is not False and len(equivs)>0:
    nToCopy = min(len(ligandsFound),len(equivs[0]))
if nToCopy > 0:
    ligandsToCopy = ligandsFound[0:nToCopy]
    coot.merge_molecules(ligandsToCopy,0)

coot.write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))'''
        self.connectSignal(self.cootPlugin,'finished',self.cootPlugin_finished)
        self.cootPlugin.process()

    @QtCore.Slot(dict)
    def cootPlugin_finished(self, status):
        print("\n\n1", status)
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(status)
        print("\n\n1","beyond")

        self.harvestFile(self.cootPlugin.container.outputData.XYZOUT[0], self.container.outputData.XYZOUT)
        #Substitute the composition section of REFMAC output to include new monomers
        #Perform analysis of output coordinate file composition
        if os.path.isfile(str(self.container.outputData.XYZOUT.fullPath)):
            from core.CCP4ModelData import CPdbData
            aCPdbData = CPdbData()
            aCPdbData.loadFile(self.container.outputData.XYZOUT.fullPath)
            #print 'aCPdbData',aCPdbData
            #print 'aCPdbData.chains',aCPdbData.composition.chains
            #print 'aCPdbData.monomers',aCPdbData.composition.monomers
            modelCompositionNode = None
            modelCompositionNodes = self.xmlroot.findall('.//ModelComposition')
            if len(modelCompositionNodes) > 0: modelCompositionNode = modelCompositionNodes[-1]
            else:
                refmacNodes = self.xmlroot.findall('.//REFMAC')
                if len(refmacNodes) > 0: modelCompositionNode = ET.SubElement(refmacNodes[-1],"ModelComposition")
            if modelCompositionNode is not None:
                for monomer in aCPdbData.composition.monomers:
                    monomerNode = ET.SubElement(modelCompositionNode,'Monomer',id=monomer)
        self.finishWithStatus(CPluginScript.SUCCEEDED)

    def harvestFile(self, pluginOutputItem, pipelineOutputItem):
        try:
            shutil.copyfile(str(pluginOutputItem.fullPath), str(pipelineOutputItem.fullPath))
            pipelineOutputItem.annotation = pluginOutputItem.annotation
            pipelineOutputItem.contentFlag = pluginOutputItem.contentFlag
            #print '#harvestFile',pluginOutputItem.fullPath, pluginOutputItem.contentFlag
            pipelineOutputItem.subType = pluginOutputItem.subType
        except:
            self.appendErrorReport(202,str(pluginOutputItem.fullPath)+' '+str(pipelineOutputItem.fullPath))
            self.finishWithStatus(CPluginScript.FAILED)

    def appendXML(self, changedFile, replacingElementOfType=None):
        newXML = CCP4Utils.openFileToEtree(changedFile)
        oldNodes = self.xmlroot.findall(replacingElementOfType)
        if len(oldNodes) > 0: oldNodes[0].parent().remove(oldNodes[0])
        self.xmlroot.append(newXML)
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlfile:
            CCP4Utils.writeXML(xmlfile,etree.tostring(self.xmlroot))

    def checkFinishStatus( self, statusDict,failedErrCode,outputFile = None,noFileErrCode= None):
        import os
        if len(statusDict)>0 and statusDict['finishStatus'] == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(statusDict['finishStatus'])
        try:
            assert outputFile.exists(),'Entity provided is not CDataFile or does not exist'
        except:
            self.appendErrorReport(noFileErrCode,'Expected file: '+str(outputFile))
            self.finishWithStatus(CPluginScript.FAILED)

    def finishWithStatus(self, status=CPluginScript.SUCCEEDED):
        self.flushXML()
        self.reportStatus(status)

    def flushXML(self):
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot))
