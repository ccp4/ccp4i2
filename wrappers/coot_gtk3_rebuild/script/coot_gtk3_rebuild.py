
import os
import glob

from lxml import etree

from ....core import CCP4PluginScript
from ....core import CCP4Utils


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

        with open(f"{self.getWorkDirectory()}/startupScript.py", "w") as outscript:
            outscript.write(f"print('{CCP4PluginScript.__file__}')\n")
            outscript.write("import sys\n")
            outscript.write("import os\n")
            outscript.write("import traceback\n")
            outscript.write(
                f"CCP4I2_ROOT='{os.path.dirname(os.path.abspath(os.path.dirname(CCP4PluginScript.__file__)))}'\n")
            outscript.write(f"sys.path.append(CCP4I2_ROOT)\n")
            outscript.write(
                f"sys.path.append(os.path.join(CCP4I2_ROOT,'wrappers','coot_gtk3_rebuild','script'))\n")
            outscript.write("try:\n")
            outscript.write(
                "    from wrappers.coot_gtk3_rebuild.script.CCP4i2 import CCP4i2CootGTK3Plugin\n")
            outscript.write(
                f"    a=CCP4i2CootGTK3Plugin.CCP4i2Menu(jobId='{self.jobId}')\n")
            outscript.write("except Exception as err:\n")
            outscript.write("    print(err)\n")
            outscript.write("    traceback.print_exc()\n")
        self.appendCommandLine(
            ['--script', os.path.join(self.getWorkDirectory(), 'startupScript.py')])

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
            self.xmlroot = etree.Element('coot_rebuild')
            e = etree.Element('number_output_files')
            e.text = str(self.numberOfOutputFiles())
            e = etree.Element('number_output_dicts')
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

                    ligNodes = self.xmlroot.xpath('//LIGANDS')
                    if len(ligNodes) == 0:
                        ligNode = etree.SubElement(self.xmlroot, 'LIGANDS')
                    else:
                        ligNode = ligNodes[0]
                    try:
                        annotation = 'Coot/Prodrg created geometry for'
                        for item in dictFile.fileContent.monomerList:
                            lig = etree.SubElement(ligNode, 'ligand')
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
        warningsNodes = self.xmlroot.xpath('//Warnings')
        if len(warningsNodes) == 0:
            warningsNode = etree.SubElement(self.xmlroot, 'Warnings')
        else:
            warningsNode = warningsNodes[0]
        warningNode = etree.SubElement(warningsNode, 'Warning')
        warningNode.text = text
