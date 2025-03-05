import os
from core.CCP4PluginScript import CPluginScript
from core.CCP4Modules import PROJECTSMANAGER


class coot1(CPluginScript):
    TASKMODULE = 'model_building'
    TASKTITLE = 'Coot 1'
    TASKNAME = 'coot1'
    TASKCOMMAND = '/home/paul/programs/ccp4-9/coot_py3/bin/coot'
    TASKVERSION = 0.0
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    ERROR_CODES = {}

    def makeCommandAndScript(self):
        script = ["import coot"]
        if self.container.inputData.XYZIN_LIST.isSet():
            for path in self.container.inputData.XYZIN_LIST:
                script.append(f"coot.read_coordinates('{path}')")
        if self.container.inputData.FPHIIN_LIST.isSet():
            for path in self.container.inputData.FPHIIN_LIST:
                script.append(f"coot.read_mtz('{path}', 'F', 'PHI', '', False, False)")
        if self.container.inputData.DELFPHIIN_LIST.isSet():
            for path in self.container.inputData.DELFPHIIN_LIST:
                script.append(f"coot.read_mtz('{path}', 'F', 'PHI', '', False, True)")
        if self.container.inputData.DELFPHIINANOM_LIST.isSet():
            for path in self.container.inputData.DELFPHIINANOM_LIST:
                script.append(f"imap = coot.read_mtz('{path}', 'F', 'PHI', '', False, True)")
                script.append("coot.set_map_colour(imap, 0.75, 0.9, 0.75)")
        jobDirectory = PROJECTSMANAGER().db().jobDirectory(self.jobId)
        scriptPath = os.path.join(jobDirectory, "coot_script.py")
        with open(scriptPath, "w") as stream:
            stream.write("\n".join(script))
        self.appendCommandLine(["--script", scriptPath])
        if self.container.inputData.COOTSCRIPTFILE.isSet():
            self.appendCommandLine(["--script", self.container.inputData.COOTSCRIPTFILE])
        return CPluginScript.SUCCEEDED
