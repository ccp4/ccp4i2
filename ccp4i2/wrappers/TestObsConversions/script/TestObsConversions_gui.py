"""
Martin Noble
"""

from PySide2 import QtCore

from ....qtgui.CCP4TaskWidget import CTaskWidget


class TestObsConversions_gui(CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'TestObsConversions'
    TASKVERSION = 0.0
    TASKMODULE='developer_tools'
    TASKTITLE="Test CCP4i2 observed data interconversions"
    WHATNEXT = []
    DESCRIPTION = '''Exercise CCP4i2 Observed data representation conversions'''
    
    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)
    
    def drawContents(self):
        folder = self.openFolder(folderFunction='inputData',title='Script control')
        self.createLine( ['label','Using file with representation','widget','INPUT_REPRESENTATION'])
        self.container.controlParameters.INPUT_REPRESENTATION.dataChanged.connect( self.handleModeChanged)
        self.createLine( [ 'label', 'Reflection object containing I+/I-','widget', '-browseDb', True, 'F_SIGF_AS_IPAIR'] ,toggle=['INPUT_REPRESENTATION','open',[1]])
        self.createLine( [ 'label', 'Reflection object containing F+/F-','widget', '-browseDb', True, 'F_SIGF_AS_FPAIR'] ,toggle=['INPUT_REPRESENTATION','open',[2]])
        self.createLine( [ 'label', 'Reflection object containing Imean','widget', '-browseDb', True, 'F_SIGF_AS_IMEAN'] ,toggle=['INPUT_REPRESENTATION','open',[3]])
        self.createLine( ['label','Representation of file to generate','widget','OUTPUT_REPRESENTATION'])
        self.container.controlParameters.OUTPUT_REPRESENTATION.dataChanged.connect( self.handleOutputModeChanged)
        self.handleModeChanged()
        self.handleOutputModeChanged()
    
    @QtCore.Slot()
    def handleModeChanged(self):
        if self.container.controlParameters.INPUT_REPRESENTATION == 1:
            self.container.controlParameters.OUTPUT_REPRESENTATION.setQualifiers({'enumerators':[1,2,3,4],'menuText':['I+ and I-', 'F+ and F-', 'Imean', 'Fmean']})
            self.container.inputData.F_SIGF_AS_IPAIR.setQualifiers({'allowUndefined':False})
            self.container.inputData.F_SIGF_AS_FPAIR.setQualifiers({'allowUndefined':True})
            self.container.inputData.F_SIGF_AS_IMEAN.setQualifiers({'allowUndefined':True})
        elif self.container.controlParameters.INPUT_REPRESENTATION == 2:
            self.container.controlParameters.OUTPUT_REPRESENTATION.setQualifiers({'enumerators':[2,4],'menuText':['F+ and F-', 'Fmean']})
            self.container.inputData.F_SIGF_AS_IPAIR.setQualifiers({'allowUndefined':True})
            self.container.inputData.F_SIGF_AS_FPAIR.setQualifiers({'allowUndefined':False})
            self.container.inputData.F_SIGF_AS_IMEAN.setQualifiers({'allowUndefined':True})
        elif self.container.controlParameters.INPUT_REPRESENTATION == 3:
            self.container.controlParameters.OUTPUT_REPRESENTATION.setQualifiers({'enumerators':[3,4],'menuText':['Imean', 'Fmean']})
            self.container.inputData.F_SIGF_AS_IPAIR.setQualifiers({'allowUndefined':True})
            self.container.inputData.F_SIGF_AS_FPAIR.setQualifiers({'allowUndefined':True})
            self.container.inputData.F_SIGF_AS_IMEAN.setQualifiers({'allowUndefined':False})
        self.getWidget('OUTPUT_REPRESENTATION').populateComboBox(self.container.controlParameters.OUTPUT_REPRESENTATION)
        self.updateViewFromModel()

    @QtCore.Slot()
    def handleOutputModeChanged(self):
        self.container.outputData.F_SIGF_FINAL.contentFlag = int(self.container.controlParameters.OUTPUT_REPRESENTATION)
        self.updateViewFromModel()
