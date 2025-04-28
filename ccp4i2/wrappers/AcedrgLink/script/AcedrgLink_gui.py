from ....qtgui.CCP4TaskWidget import CTaskWidget


#-------------------------------------------------------------------
class AcedrgLink_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'AcedrgLink' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    SHORTTASKTITLE='AceDRG in link generation mode'
    TASKTITLE='AceDRG in link generation mode'
    DESCRIPTION = '''AceDRG in link generation mode'''
    WHATNEXT = ['coot_rebuild']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
