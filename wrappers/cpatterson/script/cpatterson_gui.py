"""
    cpatterson_gui.py
    Copyright (C) 2014-2019 Jon Agirre & University of York
    Author: Jon Agirre

"""

from baselayer import QtGui, QtWidgets,QtCore
from qtgui.CCP4TaskWidget import CTaskWidget

class cpatterson_gui(CTaskWidget):

    TASKNAME = 'cpatterson'
    TASKVERSION = 0.1
    TASKMODULE='expt_data_utility'
    TASKTITLE='Calculate Patterson map'
    SHORTTASKTITLE='Calculate Patterson map'
    DESCRIPTION = 'Calculate Patterson map (cpatterson)'

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.setProgramHelpFile ( 'cpatterson' )

        self.openFolder ( folderFunction='inputData' )

        self.openSubFrame ( frame=[True] )
        self.createLine ( [ 'widget', 'F_SIGF' ] )
        self.closeSubFrame ( )

        self.closeFolder ( )
