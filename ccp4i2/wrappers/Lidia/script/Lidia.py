from pathlib import Path
import glob
import os
import platform
import shutil
import sys
import xml.etree.ElementTree as ET

from lxml import etree
from PySide2 import QtCore
from rdkit import Chem
import acedrg

from . import MOLSVG
from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript
from ....core.CCP4Modules import LAUNCHER, PREFERENCES


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
        LAUNCHER().launch(
            viewer=viewer,
            argList=argList,
            callBack=self.handleFinished,
            envEdit=envEdit,
            logFile=self.makeFileName('LOG')
        )
        return CPluginScript.SUCCEEDED

    @QtCore.Slot()
    def handleFinished(self):
        rootNode = ET.Element('Lidia')
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
                molNode = ET.SubElement(rootNode,'MOLDATA')
                molNode.text = etree.CDATA(molText)
            svgNode = ET.SubElement(rootNode,'SVGNode')
            
            svgNode.append(self.svgForMolFile(outputMOL))

        CCP4Utils.writeXml(rootNode, self.makeFileName('PROGRAMXML'))
        
        self.reportStatus(CPluginScript.SUCCEEDED)

    def svgForMolFile(self, molFilePath):
        try:
            mol = Chem.MolFromMolFile(molFilePath)
            Chem.SanitizeMol(mol)
            Chem.Kekulize(mol)
            return acedrg.svgFromMol(mol)
        except:
            return self.mySvgForMolFile(molFilePath)

    def mySvgForMolFile(self, molFilePath):
        mdlMolecule = MOLSVG.MDLMolecule(molFilePath)
        return mdlMolecule.svgXML(size=(300,300))


def _lidiaPath() -> str:
    if hasattr(PREFERENCES(), 'COOT_EXECUTABLE'):
        path = Path(str(PREFERENCES().COOT_EXECUTABLE))
        if path.is_file():
            return (str(path.resolve().parent / "lidia"),None)
    if lidiaPath := shutil.which('lidia'):
        return (str(Path(lidiaPath).resolve()),None)
    if sys.platform == "linux":# Seems that lidia does not run without PYTHONPATH being set on Linux
        return (str(Path(os.environ["CCP4"]).resolve() / "coot_py2/bin/lidia"),str(Path(os.environ["CCP4"]).resolve() / "coot_py2/"))
    return (str(Path(os.environ["CCP4"]).resolve() / "coot_py2/bin/lidia"),None)
