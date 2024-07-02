"""
    wrappers/mosflm/script/mosflm_gui.py
    Martin Noble
    """

from PySide6 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class CTaskpointless(CCP4TaskWidget.CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'mosflm'
    TASKVERSION = 0.0
    TASKMODULE='test'
    TASKTITLE='Integrate images - MOSFLM'
    DESCRIPTION='Use a script that you provide'
    WHATNEXT = ['aimless_pipe']
    
    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)
    
    
    def drawContents(self):
        
        self.setProgramHelpFile('mosflm')
        
        folder = self.openFolder(folderFunction='controlParameters',title='Control parameters',followFrom=False)
        
        self.createLine( [ 'widget', '-guiMode','multiLine','SCRIPT' ] )

