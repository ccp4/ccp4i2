"""
     ccp4mg_edit_nomrbump task widget
"""

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Modules

#-------------------------------------------------------------------
class Cccp4mg_edit_nomrbump(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'ccp4mg_edit_nomrbump'
  TASKVERSION = 0.1
  TASKMODULE='bioinformatics'
  TASKTITLE='Interactive selection of MR model components - CCP4mg'
  SHORTTASKTITLE='CCP4mg atom select'
  DESCRIPTION='Use CCP4mg to select components of a search model and output to i2 for MR'

  def drawContents(self):
    
    self.setProgramHelpFile('ccp4mg_edit_nomrbump')
                        
    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'advice', 'Input models:' ] )
    self.createLine( [ 'widget', 'XYZIN_LIST' ] )
 
  def isValid(self):
    import os
    #print 'Ccoot_rebuild.isValid'
    if self.getWidget('followFrom') is None: return
    followJobId = self.getWidget('followFrom').currentJobId()
    #print 'Ccoot_rebuild.isValid followFrom',followJobId
    return CTaskWidget.isValid(self)
  
