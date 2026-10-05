"""The one-model case of phaser_pipeline_phil: a search model, how many
copies, and optionally a structure already placed. The ensemble list the
pipeline runs with is built from them."""
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.phaser_pipeline_phil.script.phaser_pipeline_phil import phaser_pipeline_phil


class phaser_simple_phil(phaser_pipeline_phil):

    TASKNAME = "phaser_simple_phil"

    def validity(self):
        self.createEnsembleElements()
        return super().validity()

    def process(self):
        self.createEnsembleElements()
        return super().process()

    def runTimeValidity(self):
        # XYZIN_FIXED must have been placed in THIS crystal: its cell and
        # point group are checked against the data's (sameCrystalAs in the
        # def.xml), so a homologue straight from mrparse or a database, still
        # in its own crystal's frame, is refused. Haiku once gave the cyclin
        # chain of 6p8e (cell 62.4 67.5 187.3) as fixed for data of cell
        # 57.8 64.7 186.1; Phaser held it there and placed CDK4 against it
        # (LLG -892). A file left set while INPUT_FIXED is off is not used.
        error = super().runTimeValidity()
        inp = self.container.inputData
        name = f"{self.TASKNAME}.container.inputData.XYZIN_FIXED"
        if not bool(inp.INPUT_FIXED):
            error._errors = [r for r in error._errors if str(r.get("name", "")) != name]
        elif inp.XYZIN_FIXED.isSet():
            # The cell check skips a file with no cell; for a structure said
            # to be placed, having none is the finding (a predicted model's
            # CRYST1 is 1 A in P 1, which reads as no cell)
            try:
                content = inp.XYZIN_FIXED.getFileContent()
                cell = getattr(content, "cell", None) if content is not None else None
            except Exception:
                cell = None
            if cell is None:
                error.append(klass=self.TASKNAME, code=221, name=name,
                             details=("The structure given as already placed has no crystal cell "
                                      "(a predicted model, or one cut out of its entry), so it has "
                                      "not been placed in this crystal. Search for it instead: as "
                                      "the search model here, or with phaser_pipeline_phil."),
                             severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    def checkInputData(self):
        invalid = super().checkInputData()
        if not self.container.inputData.INPUT_FIXED and "XYZIN_FIXED" in invalid:
            invalid.remove("XYZIN_FIXED")
        return invalid

    def createEnsembleElements(self):
        inp = self.container.inputData
        ensembles = inp.ENSEMBLES
        ensembles.clear()
        inp.FIXENSEMBLES.clear()
        if not inp.XYZIN.isSet():
            return
        ensembles.append(ensembles.makeItem())
        search = ensembles[-1]
        search.label.set("SearchModel")
        search.number.set(int(inp.NCOPIES) if inp.NCOPIES.isSet() else 1)
        search.use.set(True)
        item = search.pdbItemList.makeItem()
        search.pdbItemList.append(item)
        item.structure.set(inp.XYZIN)
        if str(inp.ID_RMS) == "RMS":
            item.rms_to_target.set(float(inp.SEARCHRMS))
        else:
            item.identity_to_target.set(float(inp.SEARCHSEQUENCEIDENTITY))
        if inp.INPUT_FIXED.isSet() and bool(inp.INPUT_FIXED) and inp.XYZIN_FIXED.isSet():
            ensembles.append(ensembles.makeItem())
            fixed = ensembles[-1]
            fixed.label.set("KnownStructure")
            fixed.number.set(0)
            fixed.use.set(False)
            item = fixed.pdbItemList.makeItem()
            fixed.pdbItemList.append(item)
            item.structure.set(inp.XYZIN_FIXED)
            item.identity_to_target.set(0.9)
            inp.FIXENSEMBLES.append("KnownStructure")
