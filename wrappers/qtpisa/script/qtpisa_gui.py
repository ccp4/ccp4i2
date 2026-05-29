"""
     qtpisa task widget
"""

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Modules

#-------------------------------------------------------------------
class Cqtpisa(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'qtpisa'
  TASKVERSION = 0.1
  TASKMODULE='validation'
  TASKTITLE='Interface and quaternary structure analysis - PISA'
  SHORTTASKTITLE='PISA'
  DESCRIPTION='Interface and assembly analysis (qtpisa)'

  def drawContents(self):
    
    self.setProgramHelpFile('qtpisa')
                        
    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'widget', 'XYZIN' ] )
 
  def isValid(self):
    import os
    if self.getWidget('followFrom') is None: return
    followJobId = self.getWidget('followFrom').currentJobId()
    return CTaskWidget.isValid(self)
  
