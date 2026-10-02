from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile


class comit(CPluginScript):
    TASKNAME = "comit"
    TASKCOMMAND = "comit"

    def processInputFiles(self):
        # comit reads amplitudes (F_SIGF_F, F_SIGF_SIGF): convert intensities.
        # (Data from aimless are intensities; the job died in clipper with
        # "Missing column ... F_SIGF_F".)
        self.makeHklinGemmi([{"name": "F_SIGF", "target_contentFlag": CObsDataFile.CONTENT_FLAG_FMEAN},
                             "F_PHI_IN"])

    def makeCommandAndScript(self):
        params = self.container.controlParameters
        self.appendCommandLine(["-mtzin", self.workDirectory / "hklin.mtz"])
        self.appendCommandLine(["-mtzout", self.workDirectory / "hklout.mtz"])
        self.appendCommandLine(["-colin-fo", "F_SIGF_F,F_SIGF_SIGF"])
        self.appendCommandLine(["-colin-fc", "F_PHI_IN_F,F_PHI_IN_PHI"])
        self.appendCommandLine(["-colout", "i2"])
        self.appendCommandLine(["-nomit", params.NOMIT])
        self.appendCommandLine(["-pad-radius", params.PAD_RADIUS])

    def processOutputFiles(self):
        self.container.outputData.F_PHI_OUT.annotation = "Composite omit map (comit)"
        return self.splitHklout(["F_PHI_OUT"], ["i2.F_phi.F,i2.F_phi.phi"])
