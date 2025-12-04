"""
    wrappers/ProvideTLS/script/ProvideTLS_gui.py
    Martin Noble
    """

from baselayer import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets


class CProveideTLS(CCP4TaskWidget.CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'ProvideTLS'
    TASKVERSION = 0.0
    TASKMODULE='refinement'
    TASKTITLE = 'Import and/or edit TLS set definitions'     # A short title for gui menu
    SHORTTASKTITLE='Import TLS'
    DESCRIPTION = '''Enter TLS information to be used later in the project'''
    WHATNEXT = []
    
    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)
        self.becauseTLSSet = False
        self.becausePDBSet = False
    
    def drawContents(self):
        
        self.setProgramHelpFile('ProvideTLS')
        
        folder = self.openFolder(folderFunction='inputData',title='Optional objects from which to start definition')
        
        self.createLine( [ 'advice',' '] )
        self.createLine( [ 'advice',' '] )
        self.createLine( [ 'advice','Input objects from which to infer (coordinates) or copy (TLS) sets'] )
        self.createLine( [ 'widget', '-browseDb', True, 'XYZIN', 'tip', 'input coordinate set' ] )
        self.container.inputData.XYZIN.dataChanged.connect(self.handleSelectXyzin)
                         
        self.createLine( [ 'widget', '-browseDb', True, 'TLSIN', 'tip', 'starting TLS set' ] )
        self.container.inputData.TLSIN.dataChanged.connect(self.handleSelectTlsin)
                     
        self.createLine( [ 'widget', '-guiMode','multiLine','TLSTEXT' ] )

    @QtCore.Slot()
    def handleSelectXyzin(self):
        if self.becausePDBSet:
            self.becausePDBSet = False
        elif self.container.inputData.XYZIN.isSet():
            self.becausePDBSet = True
            ranges = self.rangesFromPDB( self.container.inputData.XYZIN.fullPath.__str__() )
            self.container.controlParameters.TLSTEXT = ''
            for range in ranges:
                formattedHeader = 'TLS For '+range['Type']+' chain ' + range['Chain'] + '\n'
                self.container.controlParameters.TLSTEXT += formattedHeader
                formattedRange = "RANGE '%s%4d.' '%s%4d.' \n\n" % (range['Chain'], range['firstResidue'], range['Chain'], range['lastResidue'])
                self.container.controlParameters.TLSTEXT += formattedRange
            self.container.inputData.TLSIN.unSet()

                                                         
    @QtCore.Slot()
    def handleSelectTlsin(self):
        import os
        if self.becauseTLSSet:
            self.becauseTLSSet = False
        elif self.container.inputData.TLSIN.isSet():
            self.becauseTLSSet = True
            if os.path.isfile(self.container.inputData.TLSIN.fullPath.__str__()):
                with open(self.container.inputData.TLSIN.fullPath.__str__(),'r') as myFile:
                    content = myFile.read()
                    self.container.controlParameters.TLSTEXT = content
            self.container.inputData.XYZIN.unSet()

    def rangesFromPDB(self, filePath):
        ranges = []
        import os
        if os.path.isfile(filePath):
            import mmut
            from core.CCP4ModelData import CPdbData
            aCPdbData = CPdbData()
            aCPdbData.loadFile(filePath)
            for chain in aCPdbData.composition.peptides:
                ranges.append({'Type':'Peptide', 'Chain':chain, 'firstResidue':0, 'lastResidue':999})
            for chain in aCPdbData.composition.nucleics:
                ranges.append({'Type':'Nucleic', 'Chain':chain, 'firstResidue':0, 'lastResidue':999})
        return ranges


