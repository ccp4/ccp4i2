from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile


class pairef(CPluginScript):
    TASKNAME = "pairef"
    TASKCOMMAND = "pairef"
    PERFORMANCECLASS = "CPairefPerformance"

    def processInputFiles(self):
        FMEAN = CObsDataFile.CONTENT_FLAG_FMEAN
        IMEAN = CObsDataFile.CONTENT_FLAG_IMEAN
        cols = [{"name": "F_SIGF", "target_contentFlag": FMEAN}, "FREERFLAG"]
        if self.container.inputData.F_SIGF.canConvertTo(IMEAN):
            cols += [{"name": "F_SIGF", "target_contentFlag": IMEAN}]
        self.makeHklinGemmi(cols)

    def makeCommandAndScript(self):
        inp = self.container.inputData
        par = self.container.inputParameters
        self.appendCommandLine(["--XYZIN", inp.XYZIN])
        self.appendCommandLine(["--HKLIN", self.workDirectory / "hklin.mtz"])
        if par.FIXED_TLS:
            self.appendCommandLine(["--TLSIN", inp.TLSIN])
            self.appendCommandLine(["--tls-ncyc", par.TLSCYC])
        if inp.UNMERGED.isSet():
            self.appendCommandLine(["--unmerged", inp.UNMERGED])
        if par.INIRES > 0.0:
            self.appendCommandLine(["-i", par.INIRES])
        if par.USE_PREREF:
            self.appendCommandLine(["--prerefinement-ncyc", par.NPRECYCLES])
            if par.USE_SHAKE:
                self.appendCommandLine(["--prerefinement-shake-sites", par.SHAKE])
            if par.RESETBFAC:
                self.appendCommandLine(["--prerefinement-reset-bfactor"])
        if par.SH_TYPE == "manual":
            self.appendCommandLine(["-r", par.MANSHELL])
        elif par.SH_TYPE == "semi":
            self.appendCommandLine(["-n", par.NSHELL])
            self.appendCommandLine(["-s", par.WSHELL])
        if not par.AUTO_WGT:
            self.appendCommandLine(["-w", par.WGT_TRM])
        if par.COMPLETE:
            self.appendCommandLine("--complete")
        self.appendCommandLine(["--ncyc", par.NCYCLES])
        if inp.DICT.isSet():
            self.appendCommandLine(["--libin", inp.DICT])
        if par.REFMAC_KEYWORD_FILE.isSet():
            self.appendCommandLine(["--comfile", par.REFMAC_KEYWORD_FILE])

    def processOutputFiles(self):
        path = self.workDirectory / "pairef_project" / "PAIREF_cutoff.txt"
        with path.open(encoding="utf-8") as f:
            self.container.outputData.PERFORMANCEINDICATOR.cutoff = f.read().strip()
