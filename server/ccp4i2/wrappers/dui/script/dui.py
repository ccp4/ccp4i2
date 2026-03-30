import shutil
from pathlib import Path
import gemmi

from ccp4i2.core.CCP4PluginScript import CPluginScript


class dui(CPluginScript):
    TASKNAME = 'dui'
    TASKCOMMAND = 'dui2'

    def processInputFiles(self):
        src = Path(str(self.container.inputData.DUI2_RUN_DATA))
        dst = Path(self.getWorkDirectory(), "run_dui2_nodes", "run_data")
        if src.stem == "run_data" and src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(src, dst)
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        out = self.container.outputData
        nodes = Path(self.getWorkDirectory(), "run_dui2_nodes")
        for path in sorted(nodes.glob("**/*.mtz")):
            newName = "_".join(path.relative_to(nodes).parts)
            newPath = Path(self.getWorkDirectory(), newName)
            shutil.copy(path, newPath)
            mtz = gemmi.read_mtz_file(str(path))
            outList = out.UNMERGEDMTZ if len(mtz.batches) > 0 else out.MERGEDMTZ
            outList.append(outList.makeItem())
            outList[-1].setFullPath(str(newPath))
            outList[-1].annotation = str(newName)
        return CPluginScript.SUCCEEDED
