from ccp4i2.core.CCP4PluginScript import CPluginScript


class AcedrgLink(CPluginScript):
    TASKNAME = "AcedrgLink"
    TASKCOMMAND = "acedrg"

    def makeCommandAndScript(self):
        inp = self.container.inputData
        par = self.container.controlParameters
        self.appendCommandLine(["-L", inp.INSTRUCTION_FILE])
        self.appendCommandLine(["-o", self.workDirectory / str(inp.LINK_ID)])
        if par.EXTRA_ACEDRG_KEYWORDS.isSet():
            for line in str(par.EXTRA_ACEDRG_KEYWORDS).splitlines():
                line = line.strip()
                if len(line) > 0 and line[0] != "#":
                    self.appendCommandLine(line)

    def processOutputFiles(self):
        print("AceDRG in link mode - processing output")
        inp = self.container.inputData
        out = self.container.outputData
        out.CIF_OUT.fullPath = self.workDirectory / f"{inp.LINK_ID}_link.cif"
        out.CIF_OUT.annotation = f"Link dictionary: {inp.ANNOTATION or inp.LINK_ID}"
        unl_path = self.workDirectory / f"{inp.LINK_ID}_TMP" / "UNL_for_link"
        out.UNL_PDB.fullPath = unl_path.with_suffix(".pdb")
        out.UNL_CIF.fullPath = unl_path.with_suffix(".cif")
