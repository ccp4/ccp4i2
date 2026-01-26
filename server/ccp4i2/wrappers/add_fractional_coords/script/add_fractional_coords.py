import inspect
import os

from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript


class add_fractional_coords(CPluginScript):
    TASKMODULE = "model_data_utility"
    TASKTITLE = "Add Fractional Coordinates"
    TASKNAME = "add_fractional_coords"
    TASKCOMMAND = "ccp4-python"
    MAINTAINER = "paul.bond@york.ac.uk"

    def __init__(self, *args, **kws):
        super(add_fractional_coords, self).__init__(*args, **kws)
        self.xyzin_path = os.path.join(self.getWorkDirectory(), "xyzin.pdb")
        self.xyzout_path = os.path.join(self.getWorkDirectory(), "xyzout.cif")

    def processInputFiles(self):
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.xyzin_path)
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
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
