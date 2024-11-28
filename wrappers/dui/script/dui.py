import time
import os
import re
import glob
import shutil

from core.CCP4PluginScript import CPluginScript

class dui(CPluginScript):
    TASKMODULE = 'data_processing'
    TASKTITLE = 'Integrate images - DIALS'
    DESCRIPTION = 'Launch DIALS User Interface (DUI2) and capture output'
    TASKNAME = 'dui'
    TASKCOMMAND = 'dui2'
    TASKVERSION= 0.1
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'kyle.stevenson@stfc.ac.uk'

    def __init__(self, *args, **kwargs):
        self.stime = 0.0
        self.useDialsDir = None
        self.job_dloc = None
        self.isNewRun = True
        self.bkpfileLoc = None
        CPluginScript.__init__(self, *args, **kwargs)

    def makeCommandAndScript(self):
        if str(self.container.inputData.DUI_DIR) != "":
            self.isNewRun = False
        if not self.isNewRun:
            annot = str(self.container.inputData.DUI_DIR.annotation)
            reresult = re.search(r'\(([^\)]*)', annot)
            self.job_dloc = reresult.group(1)
            self.useDialsDir = os.path.join(os.path.split(self.getWorkDirectory())[0], self.job_dloc)
            self.appendCommandLine("directory=%s"%self.useDialsDir)
        else:
            self.job_dloc = os.path.basename(os.path.normpath(self.getWorkDirectory()))
            self.useDialsDir = self.getWorkDirectory()
        self.bkpfileLoc = os.path.join(self.useDialsDir, "dui_files", "bkp.pickle") # changed from dials_files
        return CPluginScript.SUCCEEDED

    def process(self):
        self.stime = time.time()
        CPluginScript.process(self)

    def processOutputFiles(self):
        # Carry forward the bkp file - (added to carry forward loc. of original dui run). nb. take care with annot.
        if os.path.isfile(self.bkpfileLoc):
            shutil.copy(self.bkpfileLoc, self.getWorkDirectory())
            outputDIR = self.container.outputData.DUI_ODIR
            outputDIR.append(outputDIR.makeItem())
            outputDIR[-1].setFullPath(os.path.join(self.getWorkDirectory(), "bkp.pickle"))
            outputDIR[-1].annotation = "Continue from previous dials session (%s)"%self.job_dloc
        # Wipe any mtz dumped in the base directory, i2 may pick up duplicates (DUI-19.10.1)
        wipefiles = glob.glob(os.path.join(self.getWorkDirectory(), "*.mtz"))
        for wipefile in wipefiles:
            os.remove(wipefile)
        # Add the unmerged data files.
        # The file location was changed again... copy them over to where i2 expects them to be.
        newlfiles = glob.glob(os.path.join(self.useDialsDir, "dui_files", "*.mtz"))
        for afile in newlfiles:
            # Ensure file was created after this session of dui started to avoid duplicates (still safe - so keep)
            if self.stime < os.path.getmtime(afile):
                shutil.copy(afile, self.getWorkDirectory()) # be aware dui happily blots over export files.
        filelist = glob.glob(os.path.join(self.getWorkDirectory(), "*.mtz"))
        outputUNMERGED = self.container.outputData.UNMERGEDMTZ
        for afile in filelist:
            outputUNMERGED.append(outputUNMERGED.makeItem())
            outputUNMERGED[-1].setFullPath(afile)
            outputUNMERGED[-1].annotation = os.path.basename(afile)
        # Now copy the html reports to base dir as well
        newrepfiles = glob.glob(os.path.join(self.useDialsDir, "dui_files", "*.html"))
        for afile in newrepfiles:
            # Ensure file was created after this session of dui started to avoid duplicates (still safe - so keep)
            if self.stime < os.path.getmtime(afile):
                shutil.copy(afile, self.getWorkDirectory())
        return CPluginScript.SUCCEEDED

