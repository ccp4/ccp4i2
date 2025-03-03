"""
    coordinate_selector_gui.py
    Copyright (C) 2014 Newcastle University
    Author: Martin Noble
    
    """

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class coordinate_selector_gui(CTaskWidget):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'coordinate_selector'
    TASKVERSION = 0.1
    TASKMODULE = [ 'data_entry', 'model_data_utility' ]
    TASKTITLE='Import a coordinate set - optional selection of subset'
    SHORTTASKTITLE='Import coordinates'
    DESCRIPTION = '''Select (potentially complicated) subset from a coordinate set'''
    MGDISPLAYFILES = ['XYZIN','XYZOUT']
    
    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)
    
    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
        self.createLine( ['subtitle','Select a coordinate file and enter atom selection'] )
        self.createLine ( [ 'tip','Input model from which subset will be selected',
                           'widget','XYZIN' ] )
        self.getWidget('XYZIN').showAtomSelection()
        self.createLine(['tip','Specify how coordinates should be used in other tasks',
                         'label','These coordinates should be used as','widget','OVERRIDE_SUBTYPE'])
