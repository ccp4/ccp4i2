from __future__ import print_function

from core import CCP4PluginScript

import os
import glob
import shutil
from PySide2 import QtCore
from core import CCP4Modules
from core import CCP4Utils
#from lxml import etree
from xml.etree import ElementTree as ET

class coot_gtk3_rebuild(CCP4PluginScript.CPluginScript):
    # class coot_gtk3_rebuild(CInternalPlugin):

    # Where this plugin will appear on the gui
    TASKMODULE = 'model_building'
    TASKTITLE = 'Rebuild model with coot - Gtk3'     # A short title for gui menu
    # Task name - should be same as class name
    TASKNAME = 'coot_gtk3_rebuild'
    # The command to run the executable
    TASKCOMMAND = 'coot'
    TASKVERSION = 0.0                                     # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    ERROR_CODES = {200: {'description': 'Coot exited with error status'}, 201: {
        'description': 'Failed in harvest operation'}, 202: {'description': 'Failed in processOutputFiles'}}

    def makeCommandAndScript(self):
        self.dropDir = os.path.join(self.workDirectory, 'COOT_FILE_DROP')
        if not os.path.exists(self.dropDir):
            try:
                os.mkdir(self.dropDir)
            except:
                self.dropDir = self.workDirectory
                print('Could not make dropDir reset to', self.dropDir)

        self.appendCommandLine(
            [os.path.join(self.getWorkDirectory(), 'startupScript.py')])

        with open(f"{self.getWorkDirectory()}/startupScript.py", "w") as outscript:
            outscript.write(f'''
import sys
import os
import traceback
import coot_gui_api
CCP4I2_ROOT='{os.path.dirname(os.path.abspath(os.path.dirname(CCP4PluginScript.__file__)))}'
sys.path.append(CCP4I2_ROOT)
sys.path.append(os.path.join(CCP4I2_ROOT,'wrappers','coot_gtk3_rebuild','script'))
try:
    print("Hellow world 1")
    from wrappers.coot_gtk3_rebuild.script.CCP4i2 import CCP4i2CootGTK4Plugin
    print("Hellow world 2")
    a=CCP4i2CootGTK4Plugin.CCP4i2Menu(
        baseUrl='http://127.0.0.1:43434/database', 
        jobId='{self.jobId}',
        menuBar=coot_gui_api.main_toolbar(), 
        app=coot_gui_api.application()
    )
    print("Hellow world 3")
except Exception as err:
    print(err)
    traceback.print_exc()
\n''')
            

    def processOutputFiles(self):
        try:
            # First up import PDB files that have been output
            globPath = os.path.normpath(
                os.path.join(self.dropDir, 'output*.pdb'))
            outList = glob.glob(globPath)

            xyzoutList = self.container.outputData.XYZOUT
            for iFile, outputPDB in enumerate(outList):
                fpath, fname = os.path.split(outputPDB)
                iFile = int(fname[6:-4])
                os.rename(outputPDB, xyzoutList[iFile].fullPath.__str__())
                xyzoutList[iFile].annotation = "Coot output file number" + \
                    str(iFile)
                xyzoutList[iFile].subType = 1
            # Here truncate the xyzoutList back to the numberof files that werew actually found
            xyzoutList = xyzoutList[0:len(outList)]

            # Ligand builder places output cifs in the coot-cif directory as prorg-out.cif
            # 'prodrgify this residue' places output cifs in the coot-cif directory as prodrg-???.cif
            # pyrogen create "TLC"_pyrogen.cif
            cifOutList = glob.glob(os.path.normpath(
                os.path.join(self.dropDir, 'coot-ccp4', 'prodrg-*.cif')))
            cifOutList += glob.glob(os.path.normpath(os.path.join(
                self.workDirectory, 'coot-ccp4', 'prodrg-*.cif')))
            cifOutList += glob.glob(os.path.normpath(
                os.path.join(self.workDirectory, '*pyrogen.cif')))
            cifOutList += glob.glob(os.path.normpath(
                os.path.join(self.workDirectory, 'acedrg-*.cif')))

            dictoutList = self.container.outputData.DICTOUT
            for iFile, outputCIF in enumerate(cifOutList):
                fpath, fname = os.path.split(outputCIF)
                os.rename(outputCIF, dictoutList[iFile].fullPath.__str__())
                if 'acedrg' in fname:
                    annotation = 'Coot/Acedrg created geometry for ligand'
                elif 'pyrogen' in fname:
                    annotation = 'Coot/Pyrogen created geometry for ligand'
                elif 'prodrg' in fname:
                    annotation = 'Coot/Prodrg created geometry for ligand'
                else:
                    annotation = 'Coot/Prodrg created geometry for ligand'
                dictoutList[iFile].annotation = annotation
            # Here truncate the dictoutList back to the numberof files that werew actually found
            dictoutList = dictoutList[0:len(cifOutList)]

            # Create a trivial xml output file
            self.xmlroot = ET.Element('coot_rebuild')
            e = ET.Element('number_output_files')
            e.text = str(self.numberOfOutputFiles())
            e = ET.Element('number_output_dicts')
            e.text = str(len(dictoutList))
            self.xmlroot.append(e)

            # Separate out here activity to attempt merge into project dictionary....this seems flakey,
            # but is needed for ongoing work, so I ammaking it give an report a warning in case of failure, rather than
            # offer the sad face ofdoom
            try:
                for dictFile in dictoutList:
                    try:
                        self.mergeDictToProjectLib(fileName=dictFile.__str__())
                    except:
                        self.addReportWarning(
                            'mergeDictToProjectLib raised exception: Does not compromise output Dict')

                    ligNodes = self.xmlroot.findall('.//LIGANDS')
                    if len(ligNodes) == 0:
                        ligNode = ET.SubElement(self.xmlroot, 'LIGANDS')
                    else:
                        ligNode = ligNodes[0]
                    try:
                        annotation = 'Coot/Prodrg created geometry for'
                        for item in dictFile.fileContent.monomerList:
                            lig = ET.SubElement(ligNode, 'ligand')
                            lig.text = str(item.three_letter_code)
                            annotation += (' ' + str(item.three_letter_code))
                        dictFile.annotation = annotation
                    except:
                        self.addReportWarning(
                            'fileContent.monomerList raised exception: Does not compromise output Dict')
            except:
                self.addReportWarning(
                    'failed elsewhere in merging/analysing dicts: Does not compromise output Dict')
        except:
            self.appendErrorReport(202, 'Data harvesting failed')

        CCP4Utils.saveEtreeToFile(
            self.xmlroot, self.makeFileName('PROGRAMXML'))
        if (len(outList) + len(cifOutList)) > 0:
            return CCP4PluginScript.CPluginScript.SUCCEEDED
        else:
            return CCP4PluginScript.CPluginScript.MARK_TO_DELETE

    def addReportWarning(self, text):
        warningsNode = None
        warningsNodes = self.xmlroot.findall('.//Warnings')
        if len(warningsNodes) == 0:
            warningsNode = ET.SubElement(self.xmlroot, 'Warnings')
        else:
            warningsNode = warningsNodes[0]
        warningNode = ET.SubElement(warningsNode, 'Warning')
        warningNode.text = text
