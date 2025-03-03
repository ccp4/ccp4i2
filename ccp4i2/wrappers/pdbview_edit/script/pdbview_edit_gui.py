"""
     pdbview_edit task widget
"""

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Modules

#-------------------------------------------------------------------
class Cpdbview_edit(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'pdbview_edit'
  TASKVERSION = 0.1
  TASKMODULE='model_data_utility'
  TASKTITLE='Edit PDB/CIF files by hand'
  SHORTTASKTITLE='Edit PDB/CIF'
  DESCRIPTION='Edit PDB/CIF files by hand with the PdbView program'

  def drawContents(self):
    
    self.setProgramHelpFile('pdbview_edit')
                        
    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'widget', 'XYZIN_LIST' ] )
 
  def isValid(self):
    import os
    #print 'Ccoot_rebuild.isValid'
    if self.getWidget('followFrom') is None: return
    followJobId = self.getWidget('followFrom').currentJobId()
    #print 'Ccoot_rebuild.isValid followFrom',followJobId
    return CTaskWidget.isValid(self)
  
