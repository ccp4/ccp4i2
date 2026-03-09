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
  TASKVERSION = 0.1
  TASKMODULE='model_building'
  TASKTITLE='Predict regions of nucleic acids with NucleoFind'
  DESCRIPTION='Use NucleoFind to predict the regions of nucleic acid phosphates, sugars and base in an electron density map or Coulomb potential map.'
  SHORTTASKTITLE='NucleoFind'

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('nucleofind')
    folder = self.openFolder(folderFunction='inputData',title='Input Data')

    self.createLine( [ 'subtitle', 'Input map', ""] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'widget', 'FWT_PHWT' ] )
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Options', ""] )
    self.openSubFrame( frame=[True] )
    self.createLine( [ "label", "Use", "widget", "THREADS", "label", "CPU threads"])
    self.closeSubFrame()
    
    
    self.createLine( [ 'subtitle', 'Advanced Options', ""] )
    self.openSubFrame( frame=[True] )
    self.createLine( [ "label", "Truncate reflections up to", "widget", "RESOLUTION", "label", "Å resolution"])
    self.createLine( [ "widget", "SYMMETRY", "label", "Predict over entire unit cell"])
    self.createLine( [ "label", "Overlap predicted boxes by ", "widget", "OVERLAP", "label", "grid points. Lowering this may increase prediction accuracy but will increase runtime."])
    # self.createLine( [ "label", "Use GPU Acceleration", "widget", "GPU"])
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Uncertainty Options', ""] )
    self.openSubFrame( frame=[True] )
    # self.createLine( [ "label", "Uncertainity can be measured for a NucleoFind prediction either by looking at the raw probabilistic output from the deep learning model, or by computing a variance over multilpe predictions of each grid point."])
    self.createLine( [ "widget", "RAW", "label", "Use raw values",])
    self.createLine( [ "widget", "VARIANCE", "label", "Calculate and output variance maps"])
    
    self.closeSubFrame()
