#=======================================================================================
#
#    Gui Class for table one generation (KJS, 15th July 2019)
#
#=======================================================================================

from qtgui import CCP4TaskWidget

class CTasktableone(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'tableone'
    TASKVERSION = 0.0
    TASKMODULE = 'test'
    TASKTITLE = 'Generate Table One'
    SHORTTASKTITLE = 'Table One'
    DESCRIPTION = 'Generate Table One for publications.'


    def __init__(self, parent):
        CCP4TaskWidget.CTaskWidget.__init__(self, parent)


    def drawContents(self):
        self.setProgramHelpFile('tableone')
        folder = self.openFolder(folderFunction='inputData', title='Input Data and Run Parameters')
        self.createLine(['subtitle', 'Select input data', 'Observed structure factors and protein model details are required'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'UNMERGED'])
        self.createLine(['widget', 'F_SIGF'])
        self.createLine(['widget', 'FREERFLAG'])
        self.createLine(['widget', 'XYZIN'])
        self.closeSubFrame()
        return
