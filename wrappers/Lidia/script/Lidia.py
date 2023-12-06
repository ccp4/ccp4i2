from __future__ import print_function


from core.CCP4PluginScript import CPluginScript
from PySide2 import QtCore
import os,glob,re,time,sys
from core import CCP4XtalData
from lxml import etree
import math
from core import CCP4Modules
from core import CCP4Utils
import platform

class lidia(CPluginScript):
    TASKMODULE = 'wrappers'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Lidia'     # A short title for gui menu
    DESCRIPTION = 'Sketch a ligand'
    TASKNAME = 'Lidia'                                  # Task name - should be same as class name
    TASKCOMMAND = 'lidia'                                     # The command to run the executable
    if platform.system() == 'Windows': TASKCOMMAND = 'lidia.bat'
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    #WHATNEXT = ['coot_rebuild','parrot','buccaneer_build_refine_mr']
    RUNEXTERNALPROCESS=False
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    ERROR_CODES = {200 : {'description' : 'Failed to add item to mol list'},
                   201 : {'description' : 'Failed to setFullPath'},}
    
    def startProcess(self, command, **kw):
        cootExeDir = None
        if hasattr(CCP4Modules.PREFERENCES(),'COOT_EXECUTABLE'):
            if os.path.isfile(str(CCP4Modules.PREFERENCES().COOT_EXECUTABLE)):
                cootExeDir = str(CCP4Modules.PREFERENCES().COOT_EXECUTABLE)
        if cootExeDir is None:
            cootExeDir = CCP4Utils.which('coot')
        cootDir = os.path.normpath(os.path.dirname(os.path.dirname(cootExeDir)))
        envEdit = [['COOT_PREFIX',cootDir]]
        COOT_DATA_DIR =os.path.normpath(os.path.join(cootDir,'share','coot'))
        envEdit.append(['COOT_DATA_DIR',COOT_DATA_DIR])
        envEdit.append(['PWD',os.path.normpath(self.getWorkDirectory())])
        if sys.platform.startswith('linux'):
            envEdit.append(['PATH',os.path.join(os.environ["CCP4"],"libexec")])
        argList = []
        if self.container.inputData.MOLIN.isSet():
            argList.append(self.container.inputData.MOLIN.__str__())
### quick fix for 8.0.006, lidia from external coot will not work
        envEdit = [['PWD',os.path.normpath(self.getWorkDirectory())]]
        CCP4Modules.LAUNCHER().launch(viewer='lidia', argList = argList, callBack = self.handleFinished, envEdit=envEdit,logFile = self.makeFileName('LOG'))
        return CPluginScript.SUCCEEDED

    @QtCore.Slot()
    def handleFinished(self):
        rootNode = etree.Element('Lidia')
        #This is looking forward to a position where more things might work
        globPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'*.mdl'))
        outList = glob.glob(globPath)
        moloutList = self.container.outputData.MOLOUT_LIST
        for outputMOL in outList:
            fpath,fname = os.path.split(outputMOL)
            try:
                moloutList.append(moloutList.makeItem())
            except:
                self.appendErrorReport(200)
                self.reportStatus(CPluginScript.FAILED)
            try:
                moloutList[-1].setFullPath(outputMOL)
            except:
                self.appendErrorReport(201)
                self.reportStatus(CPluginScript.FAILED)
            moloutList[-1].annotation = "Lidia output file "+fname
            moloutList[-1].subType = 1
            
            with open (os.path.normpath(outputMOL),'r') as molFile:
                molText = molFile.read()
                molNode = etree.SubElement(rootNode,'MOLDATA')
                molNode.text = etree.CDATA(molText)
            svgNode = etree.SubElement(rootNode,'SVGNode')
            
            svgNode.append(self.svgForMolFile(outputMOL))
    
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            q = etree.tostring(rootNode,encoding='utf-8',pretty_print=True)
#           programXML.write(q.decode("utf-8"))
            CCP4Utils.writeXML(programXML,q)
        
        self.reportStatus(CPluginScript.SUCCEEDED)

    def svgForMolFile(self, molFilePath):
        try:
            from rdkit import Chem
            mol = Chem.MolFromMolFile(molFilePath)
            Chem.SanitizeMol(mol)
            Chem.Kekulize(mol)
            import acedrg
            return acedrg.svgFromMol(mol)
        except:
            return self.mySvgForMolFile(molFilePath)

    def mySvgForMolFile(self, molFilePath):
        from . import MOLSVG
        mdlMolecule = MOLSVG.MDLMolecule(molFilePath)
        return mdlMolecule.svgXML(size=(300,300))



