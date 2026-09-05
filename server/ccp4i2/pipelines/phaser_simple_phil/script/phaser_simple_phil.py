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
