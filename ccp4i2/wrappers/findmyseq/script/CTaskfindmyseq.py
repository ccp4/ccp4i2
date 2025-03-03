from qtgui import CCP4TaskWidget

class CTaskfindmyseq(CCP4TaskWidget.CTaskWidget):
    TASKNAME = 'findmyseq'
    TASKVERSION = 1.0
    TASKMODULE =['bioinformatics'] # Hmm not sure where bioinf probably
    TASKTITLE = 'Find My Sequence'
    SHORTTASKTITLE = "Find My Sequence"
    DESCRIPTION = 'Find the most likely sequence in a database for a given map and model'

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData', title='Input Data')
        self.createLine(['tip', 'Input reflections', 'widget', 'F_SIGF'])
        self.createLine(['tip', 'Input Map', 'widget', 'FPHI'])
        self.createLine(['tip', 'Model File', 'widget', 'XYZIN'])
        self.getWidget('XYZIN').showAtomSelection()
        self.createLine(['widget', '-browseDb', True, 'LSEQDB', 'tip', 'Local Sequence Db (see Find My Sequence instructions)'])
