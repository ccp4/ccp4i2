from __future__ import print_function

from core import CCP4PluginScript

import os
import glob
import sys
import re
import shutil
from dbapi import CCP4DbApi
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2 import QtWebEngineWidgets
from PySide2 import QtWidgets
from core import CCP4Modules
from core import CCP4Utils
from lxml import etree

class moorhen_rebuild(CCP4PluginScript.CPluginScript):
    # class moorhen_rebuild(CInternalPlugin):

    # Where this plugin will appear on the gui
    TASKMODULE = 'model_building'
    TASKTITLE = 'Rebuild model with moorhen'     # A short title for gui menu
    # Task name - should be same as class name
    TASKNAME = 'moorhen_rebuild'
    # The command to run the executable
    RUNEXTERNALPROCESS=True
    TASKVERSION = 0.0                                     # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    TASKCOMMAND = 'moorhen'

    ERROR_CODES = {200: {'description': 'Coot exited with error status'}, 201: {
        'description': 'Failed in harvest operation'}, 202: {'description': 'Failed in processOutputFiles'}}

    def makeCommandAndScript(self):
        portString = os.environ['CCP4I2_HTTP_PORT']
        self.appendCommandLine([portString, self._dbJobId])
        return CCP4PluginScript.CPluginScript.SUCCEEDED

    def processInputFiles(self):

        #def startProcess(self, *args, **kwargs):
        app = QtWidgets.QApplication.instance()    

        self.xmlroot = etree.Element('moorhen_rebuild')
        self.dropDir = os.path.join(self.workDirectory,'COOT_FILE_DROP')
        if not os.path.exists(self.dropDir):
            os.mkdir(self.dropDir)

        #rv = QtGui.QDesktopServices.openUrl(url)
        
        # This was intended for using QWebEngineView
        '''
        self.win = QtWebEngineWidgets.QWebEngineView()
        self.win.resize(1024, 768)

        self.win.load(url)

        app.setQuitOnLastWindowClosed(False)
        app.lastWindowClosed.connect(self.beforeQuit)

        self.win.show()
        self.win.raise_()
        '''
        

        self.fileSystemWatcher = QtCore.QFileSystemWatcher(parent=app)
        self.fileSystemWatcher.addPath(self.getWorkDirectory())
        self.fileSystemWatcher.directoryChanged.connect(self.handleTerminate)
        
        return CCP4PluginScript.CPluginScript.SUCCEEDED
    

    @QtCore.Slot(str)
    def handleTerminate(self, filename):
        terminateFilePath = os.path.join(self.getWorkDirectory(), "TERMINATE")
        if os.path.exists(terminateFilePath):
            print('In handle terminate', filename)
            sys.stdout.flush()
            self.postProcess()
            self.reportStatus(CCP4PluginScript.CPluginScript.SUCCEEDED)    

    def processOutputFiles(self):
        # First up import PDB files that have been pushed into the drop directory
        globPath = os.path.normpath(
            os.path.join(self.dropDir, 'output*.pdb'))
        outList = glob.glob(globPath)

        xyzoutList = self.container.outputData.XYZOUT
        for iFile, outputPDB in enumerate(outList):
            fpath, fname = os.path.split(outputPDB)
            numbers = re.findall(r'\d+', fname)
            if len(numbers) == 0:
                entryNumber = 0
            else:
                entryNumber = int(numbers[-1])

            os.rename(outputPDB, xyzoutList[entryNumber].fullPath.__str__())
            xyzoutList[entryNumber].annotation = f"Coot output file number {iFile}" 
            xyzoutList[entryNumber].subType = 1
        # Here truncate the xyzoutList back to the numberof files that werew actually found
        xyzoutList = xyzoutList[0:len(outList)]

        CCP4Utils.saveEtreeToFile(
            self.xmlroot, self.makeFileName('PROGRAMXML'))
        if len(outList) > 0:
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
