from ccp4i2.core.CCP4PluginScript import CPluginScript

class dials_rlattice(CPluginScript):

    TASKMODULE = 'data_processing'
    TASKTITLE = 'Dials_rLattice'
    TASKNAME = 'dials_rlattice'
    TASKCOMMAND = 'dials.reciprocal_lattice_viewer'
    PERFORMANCECLASS = 'CExpPhasPerformance'

    def processOutputFiles(self):
        return CPluginScript.MARK_TO_DELETE

    def makeCommandAndScript(self):
        inputJson = self.container.inputData.JSON_IN.fullPath.__str__()
        inputPickle = self.container.inputData.PICKLE_IN.fullPath.__str__()
        self.appendCommandLine(inputJson)
        self.appendCommandLine(inputPickle)
