import glob
import os
from pathlib import Path
import platform
import sys
from ccp4i2.baselayer import QtCore
from lxml import etree
from core.CCP4PluginScript import CPluginScript
from core import CCP4Modules
from core import CCP4Utils

class lidia(CPluginScript):
    TASKMODULE = 'wrappers'  # Where this plugin will appear on the gui
    TASKTITLE = 'Lidia'  # A short title for gui menu
    DESCRIPTION = 'Sketch a ligand'
    TASKNAME = 'Lidia'  # Task name - should be same as class name
    TASKCOMMAND = 'lidia.bat' if platform.system() == "Windows" else 'lidia'  # The command to run the executable
    TASKVERSION = 0.0  # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS = False
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    ERROR_CODES = {200 : {'description' : 'Failed to add item to mol list'},
                   201 : {'description' : 'Failed to setFullPath'},}
    
    def startProcess(self, command, **kw):
        viewer = 'lidia'
        argList = []
        lidiaPath = _lidiaPath()
        if not sys.platform.startswith("win"):
            viewer = '/bin/sh'
            argList = [lidiaPath[0]]
            if self.container.inputData.MOLIN.isSet():
                argList.append(str(self.container.inputData.MOLIN))
        envEdit=[['PWD', os.path.normpath(self.getWorkDirectory())]]
        if lidiaPath[1]:
            envEdit.append(["PYTHONHOME",lidiaPath[1]])
        CCP4Modules.LAUNCHER().launch(
            viewer=viewer,
            argList=argList,
            callBack=self.handleFinished,
            envEdit=envEdit,
            logFile=self.makeFileName('LOG')
        )
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


def _lidiaPath() -> str:
    if hasattr(CCP4Modules.PREFERENCES(), 'COOT_EXECUTABLE'):
        path = Path(str(CCP4Modules.PREFERENCES().COOT_EXECUTABLE))
        if path.is_file():
            return (str(path.resolve().parent / "lidia"),None)
    if lidiaPath := CCP4Utils.which('lidia'):
        return (str(Path(lidiaPath).resolve()),None)
    if sys.platform == "linux":# Seems that lidia does not run without PYTHONPATH being set on Linux
        return (str(Path(os.environ["CCP4"]).resolve() / "coot_py2/bin/lidia"),str(Path(os.environ["CCP4"]).resolve() / "coot_py2/"))
    else:
        return (str(Path(os.environ["CCP4"]).resolve() / "coot_py2/bin/lidia"),None)
