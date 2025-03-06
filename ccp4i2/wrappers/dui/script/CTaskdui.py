
from qtgui import CCP4TaskWidget

class CTaskdui(CCP4TaskWidget.CTaskWidget):

    TASKMODULE = 'data_processing'
    TASKNAME = 'dui'
    TASKVERSION = 0.01
    SHORTTASKTITLE = 'DIALS User Interface'
    TASKTITLE = 'Integrate images with DIALS'
    DESCRIPTION = 'Launch DUI and capture output'
    WHATNEXT = ['aimless_pipe']

    def __init__(self, parent):
        CCP4TaskWidget.CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.setProgramHelpFile('dui')
        self.openFolder(folderFunction='inputData', title='Input Data and Run Parameters', followFrom=False)
        self.createLine(['subtitle', 'Start the DIALS User Interface (DUI) and capture data on output'] )
        self.createLine(['label', ''])
        self.createLine(['label', 'If you wish to continue from a previous DUI session (from this project) then please'])
        self.createLine(['label', 'select continue from a previous dials session below'])
        self.createLine(['label', 'nb. do not select a dials pick file here'])
        self.createLine(['label', ''])
        self.createLine(['widget', 'DUI_DIR' ] )
        self.createLine(['label', ''])
        self.createLine(['label', 'Click "Run" to start'])
        self.createLine(['label', 'When you have finished integration in DIALS, exit DUI'])
        self.createLine(['label', 'Then use the "Data reduction" follow-on task'])

