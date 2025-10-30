from qtgui.CCP4TaskWidget import CTaskWidget


class modelASUCheck_gui(CTaskWidget):
    TASKNAME = 'modelASUCheck'
    TASKVERSION = 0.1
    TASKMODULE = ['model_data_utility']
    SHORTTASKTITLE = 'Check model against AU contents'
    TASKTITLE = 'Check model against AU contents'
    DESCRIPTION = 'Align sequences in model to those in AU contents to validate model.'

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)

        self.createLine(['subtitle', 'Input coordinates'])
        self.openSubFrame(frame=True)
        self.createLine (['tip', 'Input model', 'widget', 'XYZIN'])
        self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()

        self.createLine(['subtitle', 'Input AU Contents'])
        self.openSubFrame(frame=True)
        self.createLine (['tip', 'Input AU Contents', 'widget', 'ASUIN'])
        self.closeSubFrame()
