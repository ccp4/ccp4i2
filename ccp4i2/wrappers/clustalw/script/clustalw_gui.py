"""
Copyright (C) 2015 Newcastle University
Author: Martin Noble
"""

from PySide2 import QtCore

from ....qtgui.CCP4TaskWidget import CTaskWidget


#-------------------------------------------------------------------
class clustalw_gui(CTaskWidget):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKMODULE = 'bioinformatics'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Align sequences - CLUSTALW'     # A short title for gui menu
    SHORTTASKTITLE='CLUSTALW'
    TASKVERSION = 0.1
    DESCRIPTION = 'Align sequences using clustalw'
    TASKNAME = 'clustalw'                                  # Task name - should be same as class name

    def drawContents(self):
        self.openFolder(folderFunction='inputData')
        self.createLine(['label','Input sequences as ','stretch','widget','SEQUENCELISTORALIGNMENT'])
        self.container.inputData.SEQUENCELISTORALIGNMENT.dataChanged.connect(self.SEQASChanged)

        self.openSubFrame ( frame=[True], toggle=[ 'SEQUENCELISTORALIGNMENT', 'open', ['SEQUENCELIST'] ] )
        self.createLine(['widget','-listVisible',True,'SEQIN'])
        self.closeSubFrame()
        self.openSubFrame ( frame=[True], toggle=[ 'SEQUENCELISTORALIGNMENT', 'open', ['ALIGNMENT'] ] )
        self.createLine(['widget','ALIGNMENTIN'])
        self.closeSubFrame()
        self.SEQASChanged()

    @QtCore.Slot()
    def SEQASChanged(self):
        if self.container.inputData.SEQUENCELISTORALIGNMENT.__str__() == 'ALIGNMENT':
            self.container.inputData.SEQIN.setQualifiers({'minLength':0,'allowUndefined' : True })
            self.container.inputData.ALIGNMENTIN.setQualifiers({'allowUndefined' : False })
        elif self.container.inputData.SEQUENCELISTORALIGNMENT.__str__() == 'SEQUENCELIST':
            self.container.inputData.SEQIN.setQualifiers({'minLength':1,'allowUndefined' : False})
            self.container.inputData.ALIGNMENTIN.setQualifiers({'allowUndefined' : True })
        self.getWidget('SEQIN').validate()
        self.getWidget('ALIGNMENTIN').validate()
        #print 'SEQASChanged',self.container.inputData.SEQUENCELISTORALIGNMENT,'SEQIN',self.container.inputData.SEQIN.qualifiers('allowUndefined'),'ALIGNMENTIN',self.container.inputData.ALIGNMENTIN.qualifiers('allowUndefined')
