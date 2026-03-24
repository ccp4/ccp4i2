import glob
import os
import shutil
from pathlib import Path

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

    def processInputFiles(self):
        src = Path(str(self.container.inputData.DUI2_RUN_DATA))
        dst = Path(self.getWorkDirectory(), "run_dui2_nodes", "run_data")
        if src.stem == "run_data" and src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(src, dst)
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        # Wipe any mtz dumped in the base directory, i2 may pick up duplicates (DUI-19.10.1)
        wipefiles = glob.glob(os.path.join(self.getWorkDirectory(), "*.mtz"))
        for wipefile in wipefiles:
            os.remove(wipefile)
        # Add the unmerged data files.
        # The file location was changed again... copy them over to where i2 expects them to be.
        newlfiles = glob.glob(os.path.join(self.getWorkDirectory(), "dui_files", "*.mtz"))
        for afile in newlfiles:
            # Ensure file was created after this session of dui started to avoid duplicates (still safe - so keep)
            shutil.copy(afile, self.getWorkDirectory()) # be aware dui happily blots over export files.
        filelist = glob.glob(os.path.join(self.getWorkDirectory(), "*.mtz"))
        outputUNMERGED = self.container.outputData.UNMERGEDMTZ
        for afile in filelist:
            outputUNMERGED.append(outputUNMERGED.makeItem())
            outputUNMERGED[-1].setFullPath(afile)
            outputUNMERGED[-1].annotation = os.path.basename(afile)
        # Now copy the html reports to base dir as well
        newrepfiles = glob.glob(os.path.join(self.getWorkDirectory(), "dui_files", "*.html"))
        for afile in newrepfiles:
            # Ensure file was created after this session of dui started to avoid duplicates (still safe - so keep)
            shutil.copy(afile, self.getWorkDirectory())
        return CPluginScript.SUCCEEDED
