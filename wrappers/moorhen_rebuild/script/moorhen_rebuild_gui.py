"""
     moorhen_rebuild task widget
"""

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Modules

#-------------------------------------------------------------------
class Cmoorhen_rebuild(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'moorhen_rebuild'
  TASKVERSION = 0.1
  TASKMODULE='model_building'
  TASKTITLE='Interactive model building - Moorhen'
  SHORTTASKTITLE='Moorhen rebuild'
  DESCRIPTION='Interactive building (Moorhen)'

  WHATNEXT = ['prosmart_refmac','coot_rebuild', 'buccaneer_build_refine_mr' ]

  def drawContents(self):
    
    self.setProgramHelpFile('coot_rebuild')
                        
    self.openFolder(folderFunction='inputData')
    self.createLine(['subtitle','Key bindings'])
    self.openSubFrame(frame=True)
    self.createLine( [ 'label', 'Use Bernhard and Paul Key bindings','stretch','widget', 'USEKEYBINDINGS' ] )
    self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Coordinates' ] )
    self.createLine( [ 'widget', 'XYZIN_LIST' ] )
    self.createLine( [ 'subtitle', 'Electron density maps' ] )
    self.createLine( [ 'widget', 'FPHIIN_LIST' ] )
    self.createLine( [ 'subtitle', 'Difference density maps' ] )
    self.createLine( [ 'widget', 'DELFPHIIN_LIST' ] )
    self.createLine( [ 'subtitle', 'Anamolous difference density maps' ] )
    self.createLine( [ 'widget', 'DELFPHIINANOM_LIST' ] )
    self.createLine( [ 'subtitle', 'Additional data' ] )
    self.createLine( [ 'widget', 'DICT' ] )
    self.createLine( [ 'widget', 'COOTSCRIPTFILE' ] )
