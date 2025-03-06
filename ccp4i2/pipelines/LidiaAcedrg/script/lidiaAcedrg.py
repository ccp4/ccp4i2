
from core.CCP4PluginScript import CPluginScript
from PySide2 import QtCore
import os,glob,re,time,sys
from core import CCP4XtalData
from lxml import etree
import math
from core import CCP4Modules,CCP4Utils
from core import CCP4ErrorHandling

class lidiaAcedrg(CPluginScript):
    TASKNAME = 'LidiaAcedrg'            # Task name - should be same as class name
    TASKVERSION= 0.0                    # Version of this plugin
    ASYNCHRONOUS = False
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
            self.lidiaPlugin.waitForFinished = -1
            self.lidiaPlugin.process()
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'MOL':
            result = self.doAcedrg('MOL', self.container.inputData.MOLIN)
            self.finishWithStatus(result)
        elif self.container.inputData.MOLSMILESORSKETCH.__str__() == 'SMILES':
            result = self.doAcedrg('SMILES', self.container.inputData.SMILESIN)
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
        
        acedrgPlugin = self.makePluginObject('acedrg')
        acedrgPlugin.container.inputData.MOLORSMILES = inputType
        if inputType == 'MOL':
            acedrgPlugin.container.inputData.MOLIN = inputObject
        elif inputType == 'SMILES':
            acedrgPlugin.container.inputData.SMILESIN = inputObject
        acedrgPlugin.container.inputData.TLC = self.container.inputData.TLC
        acedrgPlugin.container.inputData.NRANDOM = self.container.inputData.NRANDOM
        acedrgPlugin.doAsync = False
        result = acedrgPlugin.process()
        if result != CPluginScript.SUCCEEDED:
            self.finishWithStatus(result)
        out = self.container.outputData
        import shutil

        if self.container.inputData.CONFORMERSFROM.__str__() != "RDKIT":
            if os.path.isfile(acedrgPlugin.container.outputData.DICTOUT.fullPath.__str__()):
                out.DICTOUT_LIST.append(out.DICTOUT_LIST.makeItem())
                out.DICTOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'.dict'))
                shutil.copyfile(acedrgPlugin.container.outputData.DICTOUT.fullPath.__str__(), out.DICTOUT_LIST[-1].fullPath.__str__())
                out.DICTOUT_LIST[-1].annotation = acedrgPlugin.container.outputData.DICTOUT.annotation
            else:
                self.appendErrorReport(201, 'DICTOUT')
            
            if os.path.isfile(acedrgPlugin.container.outputData.XYZOUT.fullPath.__str__()):
                out.XYZOUT_LIST.append(out.XYZOUT_LIST.makeItem())
                out.XYZOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'.pdb'))
                shutil.copyfile(acedrgPlugin.container.outputData.XYZOUT.fullPath.__str__(), out.XYZOUT_LIST[-1].fullPath.__str__())
                out.XYZOUT_LIST[-1].annotation = acedrgPlugin.container.outputData.XYZOUT.annotation
            else:
                self.appendErrorReport(201,'XYZOUT')

        if self.container.inputData.CONFORMERSFROM.__str__() != "ACEDRG":
            if os.path.isfile(acedrgPlugin.container.outputData.DICTOUT_RDKIT.fullPath.__str__()):
                out.DICTOUT_LIST.append(out.DICTOUT_LIST.makeItem())
                out.DICTOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'_RDKIT.dict'))
                shutil.copyfile(acedrgPlugin.container.outputData.DICTOUT_RDKIT.fullPath.__str__(), out.DICTOUT_LIST[-1].fullPath.__str__())
                out.DICTOUT_LIST[-1].annotation = acedrgPlugin.container.outputData.DICTOUT_RDKIT.annotation
            else:
                self.appendErrorReport(201, 'DICTOUT_RDKIT')

            if os.path.isfile(acedrgPlugin.container.outputData.XYZOUT_RDKIT.fullPath.__str__()):
                out.XYZOUT_LIST.append(out.XYZOUT_LIST.makeItem())
                out.XYZOUT_LIST[-1].fullPath = os.path.normpath(os.path.join(self.getWorkDirectory(),self.container.inputData.TLC.__str__()+'_RDKit.pdb'))
                shutil.copyfile(acedrgPlugin.container.outputData.XYZOUT_RDKIT.fullPath.__str__(), out.XYZOUT_LIST[-1].fullPath.__str__())
                out.XYZOUT_LIST[-1].annotation = acedrgPlugin.container.outputData.XYZOUT_RDKIT.annotation
            else:
                self.appendErrorReport(201, 'XYZOUT_RDKIT')

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
