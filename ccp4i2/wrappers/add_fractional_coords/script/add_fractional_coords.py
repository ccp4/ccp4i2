"""
    add_fractional_coords.py: CCP4 GUI Project

    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.

    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
"""

import inspect
import os
from core.CCP4ModelData import CPdbDataFile
from core.CCP4PluginScript import CPluginScript


class add_fractional_coords(CPluginScript):
    TASKMODULE = "model_data_utility"
    TASKTITLE = "Add Fractional Coordinates"
    TASKNAME = "add_fractional_coords"
    TASKCOMMAND = "ccp4-python"
    TASKVERSION = 0.1
    MAINTAINER = "paul.bond@york.ac.uk"
    ERROR_CODES = {
        201: {"description": "Failed to analyse output files"},
        202: {"description": "Failed applying selection to PDB file"},
    }
    PURGESEARCHLIST = [["hklin.mtz", 0], ["log_mtzjoin.txt", 0]]

    def __init__(self, *args, **kws):
        super(add_fractional_coords, self).__init__(*args, **kws)
        self.xyzin_path = os.path.join(self.getWorkDirectory(), "xyzin.pdb")
        self.xyzout_path = os.path.join(self.getWorkDirectory(), "xyzout.cif")

    def processInputFiles(self):
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.xyzin_path)
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, **kw):
        current_path = os.path.abspath(inspect.getfile(inspect.currentframe()))
        script_path = os.path.join(os.path.dirname(current_path), "script.py")
        self.appendCommandLine([script_path, self.xyzin_path, self.xyzout_path])
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        inputData = self.container.inputData
        outputData = self.container.outputData
        outputData.XYZOUT.setFullPath(self.xyzout_path)
        outputData.XYZOUT.annotation.set("mmCIF with fractional coordinates")
        outputData.XYZOUT.subType.set(inputData.XYZIN.subType)
        outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)
        return CPluginScript.SUCCEEDED
