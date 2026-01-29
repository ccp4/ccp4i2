

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class CTaskAreaimolUI(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'areaimol'
  TASKVERSION = 0.0
  TASKMODULE='model_data_utility'
  TASKTITLE='Solvent accessible surface - Areaimol'
  SHORTTASKTITLE='AREAIMOL'
  WHATNEXT = []
  PROGRAMHELP = 'areaimol'
  DESCRIPTION='Solvent accessible surface calculation the areaimol program'
  
  def drawContents(self):

    self.setProgramHelpFile('areaimol')
    
    folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)

    self.createLine( [ 'subtitle', 'Input Structure' ] )
    self.createLine( ['tip', 'This is the structure which the area will be calculated.', 'widget', 'XYZIN'] )
    self.createLine( [ 'subtitle', 'Keywords' ] )
    self.createLine( ['tip', 'Script'] )
    self.createLine( ['widget', '-guiMode','multiLine', 'EXTRA_AREAIMOL_KEYWORDS'] )
