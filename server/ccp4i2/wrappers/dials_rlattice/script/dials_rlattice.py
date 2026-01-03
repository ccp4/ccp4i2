from ccp4i2.core.CCP4PluginScript import CPluginScript

class dials_rlattice(CPluginScript):

    TASKMODULE = 'data_processing'      # Gui menu location
    TASKTITLE = 'Dials_rLattice'             # Short title for Gui
    TASKNAME = 'dials_rlattice'              # Task name - same as class name
    TASKCOMMAND = 'dials.reciprocal_lattice_viewer'    # The command to run the executable
    TASKVERSION = 1.0                  # plugin version
    COMTEMPLATE = None                 # The program com file template
    COMTEMPLATEFILE = None             # Name of file containing com file template
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ASYNCHRONOUS = False
    MAINTAINER = 'Kyle.Stevenson@stfc.ac.uk'

    def __init__(self, *args, **kwargs):
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        CPluginScript.process(self)

    def processInputFiles(self):
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        return CPluginScript.MARK_TO_DELETE

    def makeCommandAndScript(self, container=None):
        inputJson = self.container.inputData.JSON_IN.fullPath.__str__()
        inputPickle = self.container.inputData.PICKLE_IN.fullPath.__str__()
        self.appendCommandLine(inputJson)
        self.appendCommandLine(inputPickle)

