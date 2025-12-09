from __future__ import print_function


from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.baselayer import QtCore
import os,glob,re,time,sys,shutil
from ccp4i2.core import CCP4XtalData
from lxml import etree
import math
from ccp4i2.core import CCP4Modules
from ccp4i2.core import CCP4Utils
from ccp4i2.core import CCP4ErrorHandling

class PrepareDeposit(CPluginScript):
    TASKNAME = 'PrepareDeposit'
    TASKVERSION= 0.0
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 240
    MAXNJOBS = 4
    SUBTASKS=['refmac']
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    
    ERROR_CODES = {  200 : { 'description' : 'Alignments yielded nBestPairs != 1' },
                    201 : { 'description' : 'Failed to copy files to destination directory...do you have write access ?' },}

    def process(self):
        
        self.xmlroot = etree.Element('PrepareDeposit')
        
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()
        
        #Do sequence alignments if needed
        if self.container.inputData.PROVIDESEQUENCES:
          from ccp4i2.core import CCP4ModelData
          chainMatch = CCP4ModelData.CChainMatch(self.container.inputData.XYZIN,self.container.inputData.ASUIN)
          self.xmlroot.append(chainMatch.reportXmlAlignments())
        
        self.flushXML()
        
        refmacPlugin = self.makePluginObject('refmac')
        for propertyName in ['F_SIGF','FREERFLAG','XYZIN','TLSIN','DICT']:
            if hasattr(self.container.inputData,propertyName) and getattr(self.container.inputData,propertyName).isSet():
                setattr(refmacPlugin.container.inputData, propertyName, getattr(self.container.inputData,propertyName))
        
        #Here spoof forcing anomalous.  Really refmac should have this as explicit settable stuff
        if self.container.inputData.USINGIORF.__str__() == 'I':
            refmacPlugin.container.controlParameters.USE_TWIN.set(True)
            # Refmac now assumes to use Is if present if doing twinned refinement
            # refmacPlugin.container.controlParameters.TWIN_TYPE.set('I')
        elif self.container.inputData.USINGIORF.__str__() == 'FANOM':
            refmacPlugin.container.controlParameters.USEANOMALOUSFOR='OUTPUTMAPS'
        refmacPlugin.container.controlParameters.NCYCLES.set(0)
        
        refmacContentFlag=4
        if self.container.inputData.USINGIORF.__str__() == 'IANOM': refmacContentFlag=1
        if self.container.inputData.USINGIORF.__str__() == 'FANOM': refmacContentFlag=2
        if self.container.inputData.USINGIORF.__str__() == 'I': refmacContentFlag=3

        #Fold TLS into the B-factors ofthe output PDB
        if refmacPlugin.container.inputData.TLSIN.isSet() or refmacPlugin.container.inputData.get('AUTOTLS', False):
            refmacPlugin.container.controlParameters.TLSOUT_ADDU.set(True)
            if refmacPlugin.container.inputData.TLSIN.isSet():
                refmacPlugin.container.controlParameters.TLSMODE.set('FILE')
        refmacPlugin.container.controlParameters.NTLSCYCLES.set(0)
        refmacPlugin.container.controlParameters.NTLSCYCLES_AUTO.set(0)
        refmacPlugin.container.controlParameters.BFACSETUSE.set(False)
        refmacPlugin.container.controlParameters.B_REFINEMENT_MODE=self.container.inputData.B_REFINEMENT_MODE
    
        #Don't dip out just because DICT not available
        refmacPlugin.container.controlParameters.MAKE_NEW_LIGAND_EXIT.set(False)
        rv = refmacPlugin.process()
        
        if rv is not CPluginScript.SUCCEEDED:
            self.reportStatus(rv)

        refmacRootNode = CCP4Utils.openFileToEtree(refmacPlugin.makeFileName('PROGRAMXML'))
        self.xmlroot.append(refmacRootNode)
        self.flushXML()
        
        #Create the (potentially merged) file from which to generate the data that will be spat into the PDB
        pathForExperimentalDataToCifify = os.path.join(refmacPlugin.getWorkDirectory(),'hklin.mtz')
        if self.container.inputData.F_SIGF.contentFlag != refmacContentFlag:
            pathForExperimentalDataToCifify = os.path.join(self.getWorkDirectory(),'mergedForMakingCif.mtz')
            self.makeHklin0(miniMtzsIn=[['F_SIGF', int(self.container.inputData.F_SIGF.contentFlag)], ['F_SIGF', refmacContentFlag],['FREERFLAG',0]], hklin='mergedForMakingCif', ignoreErrorCodes=[])
        
        #Now convert refmac input to mmcif
        self.hklin2cifPlugin = self.makePluginObject('hklin2cif')
        hklinPath = os.path.normpath(pathForExperimentalDataToCifify)
        self.hklin2cifPlugin.container.inputData.HKLIN.setFullPath(hklinPath)
        rv = self.hklin2cifPlugin.process()
        if rv != CPluginScript.SUCCEEDED: self.reportStatus(rv)
        
        if hasattr(self.container.inputData,'ISREFMACMMCIFOUTPUT') and self.container.inputData.ISREFMACMMCIFOUTPUT.isSet and self.container.inputData.ISREFMACMMCIFOUTPUT:
            self.coordsToUse = self.container.inputData.XYZIN
        else:
            self.coordsToUse = refmacPlugin.container.outputData.XYZOUT
        #First pass of pdb_extract
        pdbExtract1Plugin = self.makePluginObject('pdb_extract_wrapper')
        pdbExtract1Plugin.container.inputData.XYZIN = self.coordsToUse
        pdbExtract1Plugin.container.outputData.ENTRYDATA.setFullPath(os.path.normpath(os.path.join(pdbExtract1Plugin.getWorkDirectory(),'data_template.text')))

        rv = pdbExtract1Plugin.process()
        if rv is not CPluginScript.SUCCEEDED:
            self.reportStatus(rv)
        
        modifiedTemplatePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'modifiedTemplate.text'))
        templatePath = pdbExtract1Plugin.container.outputData.ENTRYDATA.__str__()
        print('templatePath, modifiedTemplatePath',templatePath, modifiedTemplatePath)
        with open(templatePath,'r') as blankTemplate:
            lines = blankTemplate.readlines()
            with open(modifiedTemplatePath,'w') as modifiedTemplate:
                inMolecules = False
                for line in lines:
                    
                    if "molecule_entity_id" in line:
                        inMolecules = True
                    elif "CATEGORY 3:" in line:
                        inMolecules = False
                        iEntity = 1
                        for seqObj in self.container.inputData.ASUIN.fileContent.seqList:
                          if self.container.inputData.ASUIN.isSelected(seqObj):
                            modifiedTemplate.write('<molecule_entity_id="%d">\n'%(iEntity))
                            modifiedTemplate.write('<molecule_entity_type="polypeptide(L)" >\n')
                            modifiedTemplate.write('<molecule_one_letter_sequence="\n%s">\n'%(seqObj.sequence.__str__()))
                            #if len(chainsOfSequence) > iEntity-1:
                            #    modifiedTemplate.write('<molecule_chain_id="%s">\n'%(','.join(chainsOfSequence[iEntity-1])))
                            modifiedTemplate.write('< target_DB_id=" " > (if known) \n\n\n')
                            iEntity += 1
                    if not inMolecules: modifiedTemplate.write(line)

        """
        #Second pass of pdb_extract
        pdbExtract2Plugin = self.makePluginObject('pdb_extract_wrapper')
        pdbExtract2Plugin.container.inputData.XYZIN = refmacPlugin.container.outputData.XYZOUT
        pdbExtract2Plugin.container.inputData.ENTRYDATAIN.set(modifiedTemplatePath)

        rv = pdbExtract2Plugin.process()
        if rv is not CPluginScript.SUCCEEDED:
            self.reportStatus(rv)
        """
        
        import shutil
        if self.container.inputData.OUTPUTTYPE.__str__() == "DATABASE":
            structureCifPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'Coordinates.cif'))
            reflectionCifPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'Reflections.cif'))
            self.container.outputData.CIFREFLECTIONS.setFullPath(reflectionCifPath)
            self.container.outputData.CIFCOORDINATES.setFullPath(structureCifPath)
        else:
            structureCifPath = os.path.normpath(os.path.join(self.container.inputData.OUTPUT_DIRECTORY.__str__(),'Coordinates.cif'))
            reflectionCifPath = os.path.normpath(os.path.join(self.container.inputData.OUTPUT_DIRECTORY.__str__(),'Reflections.cif'))
        try:
            #shutil.copyfile(pdbExtract2Plugin.container.outputData.CIFFILE.__str__(), structureCifPath)
            shutil.copyfile(refmacPlugin.container.outputData.CIFFILE.__str__(), structureCifPath)
            shutil.copyfile(self.hklin2cifPlugin.container.outputData.CIFFILE.__str__(), reflectionCifPath)
        except:
            self.appendErrorReport(201)
            self.reportStatus(CPluginScript.FAILED)
        self.reportStatus(CPluginScript.SUCCEEDED)

    def flushXML(self, xml=None):
        if xml is None:
            if hasattr(self,'xmlroot'): xml=self.xmlroot
        import os
        tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
        with open(tmpFilename,'w') as tmpFile:
            CCP4Utils.writeXML(tmpFile,etree.tostring(xml, pretty_print=True))
        self.renameFile(tmpFilename, self.makeFileName('PROGRAMXML'))

