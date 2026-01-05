import os
import pathlib
import shutil
import xml.etree.ElementTree as ET

import coot_headless_api

from core.CCP4PluginScript import CPluginScript
from core.CCP4ModelData import CPdbDataFile


class coot_find_waters(CPluginScript):
    TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Find waters with Coot API'     # A short title for gui menu
    TASKNAME = 'coot_find_waters'  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS = False

    def startProcess(self, command=None, handler=None, **kw):
        outFormat = "cif" if self.container.inputData.XYZIN.isMMCIF() else "pdb"
        oldFullPath = pathlib.Path(str(self.container.outputData.XYZOUT.fullPath))
        if outFormat == "cif":
            self.container.outputData.XYZOUT.setFullPath(str(oldFullPath.with_suffix('.cif')))
            self.container.outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)

        xyzin = str(self.container.inputData.XYZIN.fullPath)
        mtzin = str(self.container.inputData.FPHIIN.fullPath)
        xyzout = str(self.container.outputData.XYZOUT.fullPath)
        threshold = self.container.controlParameters.THRESHOLD
        mindist = self.container.controlParameters.MINDIST
        maxdist = self.container.controlParameters.MAXDIST
        
        mc = coot_headless_api.molecules_container_py(True)
        mc.set_make_backups(False)
        mc.set_use_gemmi(False)
        imol = mc.read_pdb(xyzin)
        imap = mc.read_mtz(mtzin, "F", "PHI", "", False, False)
        mc.set_add_waters_sigma_cutoff(threshold)
        mc.set_add_waters_water_to_protein_distance_lim_min(mindist)
        mc.set_add_waters_water_to_protein_distance_lim_min(maxdist)
        nwaters = mc.add_waters(imol, imap)
        mc.write_coordinates(imol, xyzout)

        root = ET.Element('coot_find_waters')
        element = ET.SubElement(root, "WatersFound")
        element.text = str(nwaters)
        tree = ET.ElementTree(root)
        tree.write(self.makeFileName('PROGRAMXML'))

        shutil.rmtree("coot-backup", ignore_errors=True)

        status = CPluginScript.FAILED
        if os.path.exists(xyzout):
            status = CPluginScript.SUCCEEDED
        self.reportStatus(status)
        return status
