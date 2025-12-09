"""
     ccp4mg_general task widget
"""

from ccp4i2.baselayer import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.core import CCP4Modules

#-------------------------------------------------------------------
class Cccp4mg_general(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'ccp4mg_general'
  TASKVERSION = 0.1
  TASKMODULE='model_building'
  TASKTITLE='Molecular graphics visualization and figure creation - CCP4MG'
  SHORTTASKTITLE='CCP4MG'
  DESCRIPTION='Interactive molecular graphics: visualization, figure preparation, analysis.'

  def drawContents(self):
    
    self.setProgramHelpFile('ccp4mg_general')
                        
    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'widget', 'XYZIN_LIST' ] )
    self.createLine( [ 'advice', 'Electron density maps:' ] )
    self.createLine( [ 'widget', 'FPHIIN_LIST' ] )
    self.createLine( [ 'advice', 'Difference density maps:' ] )
    self.createLine( [ 'widget', 'DELFPHIIN_LIST' ] )
    self.createLine( [ 'advice', 'Sequences:' ] )
    self.createLine( [ 'widget', 'SEQUENCE_LIST' ] )
    self.createLine( [ 'advice', 'Ligand dictionary:' ] )
    self.createLine( [ 'widget', 'DICT' ] )
 
  def isValid(self):
    import os
    #print 'Ccoot_rebuild.isValid'
    if self.getWidget('followFrom') is None: return
    followJobId = self.getWidget('followFrom').currentJobId()
    #print 'Ccoot_rebuild.isValid followFrom',followJobId
    return CTaskWidget.isValid(self)
  
