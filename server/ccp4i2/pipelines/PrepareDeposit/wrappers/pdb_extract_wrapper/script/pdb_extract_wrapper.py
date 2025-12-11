from __future__ import print_function


from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.baselayer import QtCore
import os,glob,re,time,sys
from ccp4i2.core import CCP4XtalData
from lxml import etree
import math
from ccp4i2.core import CCP4Modules
from ccp4i2.core import CCP4Utils

class pdb_extract_wrapper(CPluginScript):
    TASKMODULE = 'wrappers'                               # Where this plugin will appear on the gui
    TASKTITLE = 'pdb_extract_wrapper'     # A short title for gui menu
    DESCRIPTION = 'Run pdb_extract'
    TASKNAME = 'pdb_extract_wrapper'                                  # Task name - should be same as class name
    TASKCOMMAND = 'pdb_extract'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS=False

    ERROR_CODES = {  200 : { 'description' : 'Failed to add item to mol list' },201 : { 'description' : 'Failed to setFullPath' },}
    
    def process(self):
        self.xmlroot = etree.Element('pdb_extract_wrapper')
        
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()
        
        PDB_EXTRACT_DIR = os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'pipelines','PrepareDeposit','wrappers','pdb_extract_wrapper'))
        envEdit = [['PDB_EXTRACT',PDB_EXTRACT_DIR]]
        envEdit.append(['PWD',os.path.normpath(self.getWorkDirectory())])
        
        if self.container.inputData.ENTRYDATAIN.isSet():
            argList = ['-r','refmac5','-ipdb',self.container.inputData.XYZIN.__str__(),'-iENT',self.container.inputData.ENTRYDATAIN.__str__()]
            print('Arglist',argList)
            print('envEdit',envEdit)
            p = CCP4Modules.LAUNCHER().launch(viewer='pdb_extract', argList = argList, callBack = self.handleFinished, envEdit=envEdit,logFile = self.makeFileName('LOG'))
            p.waitForFinished(-1)
        else:
            #MN: Note the nasty kludge of using "-SOL" to force an inpfile name to be generated to work round what looks to me to be a bug in pdb_extract v3.22 code
            argList = ['-PDB',self.container.inputData.XYZIN.__str__(),'-SOL',os.path.join(self.getWorkDirectory(),'intermediateFile.pdb')]
            print('Arglist',argList)
            print('envEdit',envEdit)
            p = CCP4Modules.LAUNCHER().launch(viewer='extract', argList = argList, callBack = self.handleFinished, envEdit=envEdit,logFile = self.makeFileName('LOG'))
            p.waitForFinished(-1)
        return CPluginScript.SUCCEEDED

    @QtCore.Slot(int,int)
    def handleFinished(self,exitCode,exitStatus):
        if not self.container.inputData.ENTRYDATAIN.isSet():
            entryPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'data_template.text'))
            self.container.outputData.ENTRYDATA.setFullPath(entryPath)
        else:
            filePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'pdb_extract.mmcif'))
            self.container.outputData.CIFFILE.setFullPath(filePath)
            print('pdb_extract.mmcif',self.container.outputData.CIFFILE.__str__(), os.path.isfile(self.container.outputData.CIFFILE.__str__()))
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot, pretty_print=True))
        
        self.reportStatus(CPluginScript.SUCCEEDED)

