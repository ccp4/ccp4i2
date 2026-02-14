from ccp4i2.core.CCP4PluginScript import CPluginScript

class dials_image(CPluginScript):

    TASKTITLE = 'Dials_image'
    TASKNAME = 'dials_image'
    TASKCOMMAND = 'dials.image_viewer'
    PERFORMANCECLASS = 'CExpPhasPerformance'

    def processOutputFiles(self):
        return CPluginScript.MARK_TO_DELETE

    def makeCommandAndScript(self):
        inputJson = self.container.inputData.JSON_IN.fullPath.__str__()
        self.appendCommandLine(inputJson)
        if self.container.inputData.PICKLE_IN.isSet():
            inputPickle = self.container.inputData.PICKLE_IN.fullPath.__str__()
            self.appendCommandLine(inputPickle)
