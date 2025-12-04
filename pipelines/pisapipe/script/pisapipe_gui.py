from qtgui.CCP4TaskWidget import CTaskWidget
from baselayer import QtCore

#-------------------------------------------------------------------
class CTaskPisapipe(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'pisapipe'
  TASKVERSION = 0.1
  TASKMODULE='test'
  TASKTITLE='Structure analysis with Pisa'
  DESCRIPTION = 'Analyse tertiary structure and interfaces of a protein'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('pisapipe')

    self.openFolder(folderFunction='inputData',followFrom=False)

    self.createLine ( [ 'widget','PDBIN' ] )
    
