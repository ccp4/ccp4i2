

from PySide6 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class CTaskPDBSetUI(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'pdbset_ui'
  TASKVERSION = 0.0
  TASKMODULE='model_data_utility'
  TASKTITLE='Scripted structure edits - Pdbset'
  SHORTTASKTITLE='PDBSET'
  WHATNEXT = []
  PROGRAMHELP = 'pdbset'
  DESCRIPTION='Structure edits with the pdbset program'
  
  def drawContents(self):

    self.setProgramHelpFile('pdbset')
    
    folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)

    self.createLine( [ 'subtitle', 'Input Structure' ] )
    self.createLine( ['tip', 'This is the structure to which edits will be applied.', 'widget', 'XYZIN'] )
    self.createLine( [ 'subtitle', 'Keywords' ] )
    self.createLine( ['tip', 'Script'] )
    self.createLine( ['widget', '-guiMode','multiLine', 'EXTRA_PDBSET_KEYWORDS'] )
