from __future__ import print_function
"""
     tasks/edstats/CTaskEdstats.py: CCP4 GUI Project
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

"""
     Jon Agirre         2014 - Started development

"""
from PySide6 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets
import multiprocessing

class CTaskEdstats(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'edstats'
  TASKVERSION = 0.1
  TASKMODULE='validation'
  TASKTITLE='Analyse agreement between model and density - EDSTATS'
  SHORTTASKTITLE='EDSTATS'
  DESCRIPTION='Calculates real-space metrics for evaluating the agreement between model and density (Edstats, cfft)'
  WHATNEXT = [ 'coot_rebuild' ]
  RANK=1
  ERROR_CODES = {  200 : { 'description' : 'Model and map cell dimensions do not match' },
                    201 : { 'description' : 'Density map and difference density map cell dimensions do not match' },}

  def taskValidity(self):
      from core import CCP4ErrorHandling
      rv = CCP4ErrorHandling.CErrorReport()
      if self.container.inputData.FPHIIN1.isSet():
          if self.container.inputData.XYZIN.isSet():
              cellsAreSame = self.container.inputData.FPHIIN1.fileContent.clipperSameCellCoor(self.container.inputData.XYZIN.fileContent)
              if not cellsAreSame['validity']:
                  modelcell = self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()
                  modeltext = ""
                  modeltext = modeltext + '\nCell:'
                  for ii in range(1,7):
                     modeltext = modeltext + '  ' + str(modelcell[ii])
                  mapcell = self.container.inputData.FPHIIN1.fileContent.cell
                  maptext = str(mapcell.a) + ' ' + str(mapcell.b) + ' ' + str(mapcell.c) + ' ' + str(mapcell.alpha) + ' ' + str(mapcell.beta) + ' ' + str(mapcell.gamma)
                  rv.append(self.__class__,200,details='Model cell dimensions:\n'+ modeltext +'\nMap cell dimensions:\n'+ maptext ,stack=False)
             
          if self.container.inputData.FPHIIN2.isSet():
              cellsAreSame = self.container.inputData.FPHIIN1.fileContent.clipperSameCell(self.container.inputData.FPHIIN2.fileContent)
              if not cellsAreSame['validity']:
                  mapcell = self.container.inputData.FPHIIN1.fileContent.cell
                  maptext = str(mapcell.a) + ' ' + str(mapcell.b) + ' ' + str(mapcell.c) + ' ' + str(mapcell.alpha) + ' ' + str(mapcell.beta) + ' ' + str(mapcell.gamma)
                  mapcell2 = self.container.inputData.FPHIIN2.fileContent.cell
                  maptext2 = str(mapcell2.a) + ' ' + str(mapcell2.b) + ' ' + str(mapcell2.c) + ' ' + str(mapcell2.alpha) + ' ' + str(mapcell2.beta) + ' ' + str(mapcell2.gamma)
                  rv.append(self.__class__,201,details='Map cell dimensions::\n'+ maptext +'\nDifference map cell dimensions:\n'+ maptext2 ,stack=False)

      return rv

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('edstats')


# the input data tab starts here

    folder = self.openFolder(folderFunction='inputData',title='Input Data' )

    self.createLine( [ 'subtitle', 'Model to analyse', 'Input your model in PDB or mmCIF format' ] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'widget', 'XYZIN'])
    self.closeSubFrame ( )

    self.createLine( [ 'subtitle', 'Map coefficients' , 'Input your 2mFo-DFc and mFo-DFc map coefficients'  ] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'advice', 'Electron density map' ] )
    self.createLine( [ 'widget', 'FPHIIN1' ] )

    self.createLine( [ 'advice', 'Difference density map' ] )
    self.createLine( [ 'widget', 'FPHIIN2' ] )
    self.closeSubFrame ( )

    self.closeFolder()

    folder = self.openFolder(folderFunction='controlParameters', title='Options' )

    self.createLine( [ 'advice', 'Specify the resolution range you would like to use' ] )
    self.createLine( [ 'label', 'Low resolution limit ', 'widget', 'RES_LOW', 'stretch', 1, 'label', 'High resolution limit ' , 'widget', 'RES_HIGH'  ] )

    self.createLine( [ 'advice', 'Map values are averaged separately for main- and side-chains' ] )
    self.createLine( [ 'label', 'Main-chain averaging is to be performed across ', 'widget', 'MAIN_AVERAGING' ] )
    self.createLine( [ 'label', 'Side-chain averaging is to be performed across ', 'widget', 'SIDE_AVERAGING' ] )

    self.createLine( [ 'widget', 'SCALING', 'label', 'Rescale Q-Q plot using ', 'widget', 'SCALING_TYPE' ] )

    self.container.inputData.FPHIIN1.dataChanged.connect( self.updateLimits )

    self.createLine( [ 'advice', 'Adjust the rejection criteria for flagging up the outliers in Coot' ] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'label', 'Accuracy metrics: RZ- < ', 'widget', 'SIGMA_RZ_MINUS', 'label', ' sigma', 'stretch', 1, 'label', 'RZ+ > ', 'widget', 'SIGMA_RZ_PLUS', 'label', ' sigma' ] )
    self.createLine( [ 'label', 'Precision metric: RO < ', 'widget', 'SIGMA_RO' ] )
    self.closeSubFrame ( )
    self.createLine( [ 'widget', 'OUTPUT_PDB_FILE', 'label', 'Output PDB file with per-atom metrics' ] )

    self.closeFolder() 
    self.updateLimits ( )

  @QtCore.Slot()
  def updateLimits ( self ):

    if self.container.inputData.FPHIIN1.isSet ( ):
        self.container.inputData.FPHIIN1.loadFile()
        self.container.inputData.RES_LOW.set ( '%.2f'%( self.container.inputData.FPHIIN1.fileContent.resolutionRange.low ) )
        self.container.inputData.RES_HIGH.set ( '%.2f'%( self.container.inputData.FPHIIN1.fileContent.resolutionRange.high ) )
        self.getWidget ( 'RES_LOW' ).validate ( )
        self.getWidget ( 'RES_HIGH' ).validate ( )
    return
