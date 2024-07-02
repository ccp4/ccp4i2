"""
     tasks/coot_fit_residues
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from PySide6 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class Ccoot_fit_residues(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'coot_fit_residues'
  TASKVERSION = 0.1
  #TASKMODULE='model_building'
  TASKTITLE='Coot fit residues'

  def drawContents(self):

      
    self.setProgramHelpFile('coot_fit_residues')
                        
    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'tip', 'input structure', 'widget', 'XYZIN' ] )
    self.createLine( [ 'tip', 'input 2mFo-DFc coefficient', 'widget', 'FPHIIN' ] )
 
