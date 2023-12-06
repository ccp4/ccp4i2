"""
    wrappers/import_mosflm/script/import_mosflm_gui.py
    Martin Noble
    """

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class CTaskImportMosflm(CCP4TaskWidget.CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'import_mosflm'
    TASKVERSION = 0.0
    TASKTITLE='Import iMosflm X-ray data'
    TASKMODULE='wrappers'
    DESCRIPTION = 'Import merged and unmerged X-ray reflections from Mosflm'
    WHATNEXT = ['aimless_pipe']
    
    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)
    
    
    def drawContents(self):
        
        self.setProgramHelpFile('import_mosflm')
        
        folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)
        
        self.createLine( ['widget', 'MOSFLMXML'] )
        self.createLine( ['widget', 'UNMERGEDMTZ'] )

