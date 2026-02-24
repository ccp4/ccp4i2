from ccp4i2.core.CCP4PluginScript import CPluginScript


class comit(CPluginScript):
    TASKNAME = "comit"
    TASKCOMMAND = "comit"

    def processInputFiles(self):
        self.makeHklinGemmi(["F_SIGF", "F_PHI_IN"])

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
        return self.splitHklout(["F_PHI_OUT"], ["i2.F_phi.F,i2.F_phi.phi"])
