"""
     tasks/privateer/privateer_gui.py: CCP4 GUI Project
     Copyright (C) 2014 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the
     license to address the requirements of UK law.

     You should have received a copy of the modified GNU Lesser General
     Public License along with this library.  If not, copies may be
     downloaded from http://www.ccp4.ac.uk/ccp4license.php

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""


from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class NucleoFindGUI(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'nucleofind'
  TASKVERSION = 0.3
  TASKMODULE='model_building'
  TASKTITLE='Predict nucleic acid positions - NucleoFind'
  DESCRIPTION=''
  SHORTTASKTITLE='NucleoFind'

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('nucleofind')


# the input data tab starts here

    folder = self.openFolder(folderFunction='inputData',title='Input Data')

    self.createLine( [ 'subtitle', 'Input map', ""] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'widget', 'FWT_PHWT' ] )
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Advanced Options', ""] )
    self.openSubFrame( frame=[True] )
    self.createLine( [ "label", "Truncate reflections up to", "widget", "RESOLUTION", "label", "Å resolution"])
    self.createLine( [ "label", "Predict over entire unit cell", "widget", "SYMMETRY"])
    self.createLine( [ "label", "Overlap predicted boxes by ", "widget", "OVERLAP", "label", "grid points"])
    self.createLine( [ "label", "Use GPU Acceleration", "widget", "GPU"])
    
    self.closeSubFrame()


  @QtCore.Slot()
  def updateRequirements ( self ) :
    if self.container.controlParameters.NEW_SUGAR :
        self.container.controlParameters.CODEIN.setQualifier ( 'allowUndefined', False )
        self.container.controlParameters.CODEIN.validate()
    else :
        self.container.controlParameters.CODEIN.setQualifier ( 'allowUndefined', True )
        self.container.controlParameters.CODEIN.validate()
