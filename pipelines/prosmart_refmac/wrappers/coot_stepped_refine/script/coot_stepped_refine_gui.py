"""
     tasks/coot_stepped_refine
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from baselayer import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class Ccoot_stepped_refine(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'coot_stepped_refine'
  TASKVERSION = 0.1
  #TASKMODULE='model_building'
  TASKTITLE='Coot stepped refinement'

  def drawContents(self):

      
    self.setProgramHelpFile('coot_stepped_refine')
                        
    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'tip', 'input structure', 'widget', 'XYZIN' ] )
    self.createLine( [ 'tip', 'input 2mFo-DFc coefficient', 'widget', 'FPHIIN' ] )

    self.openFolder(folderFunction='controlParameters',title='Options')
    self.createLine( [ 'label', 'Use Ramachandran restraints: ', 'tip', 'encourage Ramachandran compliance', 'stretch', 'widget', 'USERAMA' ] )


