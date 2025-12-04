"""
    wrappers/imosflm/script/CTaskimosflm.py
    Phil Evans
    """

from baselayer import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class CTaskiMasflm(CCP4TaskWidget.CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKMODULE = 'data_processing'                               # Where this plugin will appear on the gui
    TASKNAME = 'imosflm'
    TASKVERSION = 0.0
    TASKTITLE = 'Integrate images with Mosflm'     # A short title for gui menu
    DESCRIPTION = 'Launch iMosflm and capture output'
    WHATNEXT = ['aimless_pipe']
    
    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)
    
    
    def drawContents(self):
        
        self.setProgramHelpFile('imosflm')
        folder = self.openFolder(folderFunction='inputData',title='Input Data')

        self.createLine( ['subtitle', 'Start iMosflm and capture data on output'] )
        self.createLine( ['label',''] )
        self.createLine( ['label','Click "Run" to start'] )
        self.createLine( ['label', \
           'When you have finished integration in iMosflm, select "Session->Exit" to return here'] )
        self.createLine( ['label', 'Then use the "Data reduction" follow-on task'])

