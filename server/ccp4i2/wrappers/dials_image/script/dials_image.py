from ccp4i2.core.CCP4PluginScript import CPluginScript

class dials_image(CPluginScript):

    TASKMODULE = 'data_processing'      # Gui menu location
    TASKTITLE = 'Dials_image'             # Short title for Gui
    TASKNAME = 'dials_image'              # Task name - same as class name
    TASKCOMMAND = 'dials.image_viewer'    # The command to run the executable
    PERFORMANCECLASS = 'CExpPhasPerformance'
    MAINTAINER = 'Kyle.Stevenson@stfc.ac.uk'

    def processOutputFiles(self):
        return CPluginScript.MARK_TO_DELETE

    def makeCommandAndScript(self):
        inputJson = self.container.inputData.JSON_IN.fullPath.__str__()
        self.appendCommandLine(inputJson)
        if self.container.inputData.PICKLE_IN.isSet():
            inputPickle = self.container.inputData.PICKLE_IN.fullPath.__str__()
            self.appendCommandLine(inputPickle)
