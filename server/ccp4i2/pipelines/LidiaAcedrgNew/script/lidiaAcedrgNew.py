import os
import sys

from lxml import etree
from ccp4i2.baselayer import QtCore

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class lidiaAcedrgNew(CPluginScript):
    TASKNAME = 'LidiaAcedrgNew'            # Task name - should be same as class name
    TASKVERSION= 0.0                    # Version of this plugin
    TIMEOUT_PERIOD = 9999999.9
    WHATNEXT = ['coot_rebuild']
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Expected output file not made', 'severity':CCP4ErrorHandling.SEVERITY_WARNING },}

    def process(self):
        
        self.xmlroot = etree.Element('LidiaAcedrg')
        tlcNode = etree.SubElement(self.xmlroot,'TLC')
        tlcNode.text = 'UNL'
        if self.container.inputData.TLC.isSet():
            tlcNode.text = self.container.inputData.TLC.__str__()
        
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()
        
        if self.container.inputData.MOLSMILESORSKETCH.__str__() == 'SKETCH':
            self.lidiaPlugin = self.makePluginObject('Lidia')
            if self.container.inputData.MOLIN.isSet():
                self.lidiaPlugin.container.inputData.MOLIN = self.container.inputData.MOLIN
            self.connectSignal(self.lidiaPlugin,'finished',self.lidiaFinished)
            self.lidiaPlugin.doAsync = True
            self.lidiaPlugin.process()
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'SMILESFILE':
            result = self.doAcedrg('SMILESFILE', self.container.inputData.SMILESFILEIN)
            self.finishWithStatus(result)
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'MOL':
            result = self.doAcedrg('MOL', self.container.inputData.MOLIN)
            self.finishWithStatus(result)
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'MOL2':
            result = self.doAcedrg('MOL2', self.container.inputData.MOL2IN)
            self.finishWithStatus(result)
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'SMILES':
            result = self.doAcedrg('SMILES', self.container.inputData.SMILESIN)
            self.finishWithStatus(result)
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'DICT':
            result = self.doAcedrg('DICT', self.container.inputData.DICTIN2)
            self.finishWithStatus(result)
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'PDBMMCIF':
            result = self.doAcedrg('PDBMMCIF', self.container.inputData.PDBMMCIFIN)
            self.finishWithStatus(result)

    @QtCore.Slot(dict)
    def lidiaFinished(self, statusDict):
        status = statusDict['finishStatus']
        if status == CPluginScript.FAILED:
            self.finishWithStatus(status)
        
        lidiaRoot = CCP4Utils.openFileToEtree(self.lidiaPlugin.makeFileName('PROGRAMXML'))
        self.xmlroot.append(lidiaRoot)
        self.flushXML()
        
        out = self.container.outputData
        for molFile in self.lidiaPlugin.container.outputData.MOLOUT_LIST:
            
            out.MOLOUT_LIST.append(self.container.outputData.MOLOUT_LIST.makeItem())
            out.MOLOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'.mol'))
            import shutil
            shutil.copyfile(molFile.fullPath.__str__(),out.MOLOUT_LIST[-1].fullPath.__str__())
            out.MOLOUT_LIST[-1].annotation = 'MOL file for '+self.container.inputData.TLC.__str__()
            
            result = self.doAcedrg('MOL', molFile)
            if result != CPluginScript.SUCCEEDED: self.finishWithStatus(result)
        self.finishWithStatus(CPluginScript.SUCCEEDED)

    def doAcedrg(self, inputType, inputObject):
        
        acedrgPlugin = self.makePluginObject('acedrgNew')
        acedrgPlugin.container.inputData.MOLORSMILES = inputType
        if inputType == 'MOL':
            acedrgPlugin.container.inputData.MOLIN = inputObject
        elif inputType == 'MOL2':
            acedrgPlugin.container.inputData.MOL2IN = inputObject
        elif inputType == 'SMILESFILE':
            acedrgPlugin.container.inputData.SMILESFILEIN = inputObject
        elif inputType == 'SMILES':
            acedrgPlugin.container.inputData.SMILESIN = inputObject
        elif inputType == 'DICT':
            acedrgPlugin.container.inputData.DICTIN2 = inputObject
        elif inputType == 'PDBMMCIF':
            acedrgPlugin.container.inputData.PDBMMCIFIN = inputObject
        try:
            acedrgPlugin.container.inputData.TLC.set(self.container.inputData.TLC)
            acedrgPlugin.container.controlParameters.NOPROT.set(self.container.controlParameters.NOPROT)
            acedrgPlugin.container.controlParameters.USE_COORD.set(self.container.controlParameters.USE_COORD)
            if self.container.controlParameters.TOGGLE_NRANDOM:
                acedrgPlugin.container.inputData.NRANDOM.set(self.container.inputData.NRANDOM)
            else:
                acedrgPlugin.container.inputData.NRANDOM.set(0)
            if self.container.controlParameters.TOGGLE_METAL:
                acedrgPlugin.container.inputData.METAL_STRUCTURE = self.container.inputData.METAL_STRUCTURE
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stdout.write(str(exc_type)+'\n')
            sys.stdout.write(str(exc_value)+'\n')
            raise

        myMatchTLC = self.container.inputData.MATCHTLC
        print("My MATCHTLC",myMatchTLC.__str__())

        myATOMMATCHOPTION = self.container.inputData.ATOMMATCHOPTION
        print("My ATOMMATCHOPTION",myATOMMATCHOPTION.__str__())

        myLOCALDICT = self.container.inputData.DICTIN
        print("My ATOMMATCHOPTION",myLOCALDICT.__str__())

        acedrgPlugin.container.inputData.MATCHTLC = myMatchTLC
        acedrgPlugin.container.inputData.ATOMMATCHOPTION = myATOMMATCHOPTION
        acedrgPlugin.container.inputData.DICTIN = myLOCALDICT

        result = acedrgPlugin.process()
        if result != CPluginScript.SUCCEEDED:
            self.finishWithStatus(result)
        out = self.container.outputData
        import shutil

        if os.path.isfile(acedrgPlugin.container.outputData.DICTOUT.fullPath.__str__()):
            out.DICTOUT_LIST.append(out.DICTOUT_LIST.makeItem())
            out.DICTOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'.cif'))
            shutil.copyfile(acedrgPlugin.container.outputData.DICTOUT.fullPath.__str__(), out.DICTOUT_LIST[-1].fullPath.__str__())
            out.DICTOUT_LIST[-1].annotation = acedrgPlugin.container.outputData.DICTOUT.annotation
        else:
            self.appendErrorReport(201, 'DICTOUT')

        # No output PDB for 5-letter code monomers
        if os.path.isfile(acedrgPlugin.container.outputData.XYZOUT.fullPath.__str__()):
            out.XYZOUT_LIST.append(out.XYZOUT_LIST.makeItem())
            out.XYZOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'.pdb'))
            shutil.copyfile(acedrgPlugin.container.outputData.XYZOUT.fullPath.__str__(), out.XYZOUT_LIST[-1].fullPath.__str__())
            out.XYZOUT_LIST[-1].annotation = acedrgPlugin.container.outputData.XYZOUT.annotation

        out.MOLOUT_LIST.append(out.MOLOUT_LIST.makeItem())
        out.MOLOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'_RDKIT.mol'))
        shutil.copyfile(acedrgPlugin.container.outputData.MOLOUT.fullPath.__str__(), out.MOLOUT_LIST[-1].fullPath.__str__())
        out.MOLOUT_LIST[-1].annotation = acedrgPlugin.container.outputData.MOLOUT.annotation

        acedrgRoot = CCP4Utils.openFileToEtree(acedrgPlugin.makeFileName('PROGRAMXML'))
        self.xmlroot.append(acedrgRoot)

        return CPluginScript.SUCCEEDED

    def finishWithStatus(self, status=CPluginScript.SUCCEEDED):
        self.flushXML()
        self.reportStatus(status)

    def flushXML(self):
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot,pretty_print=True))
