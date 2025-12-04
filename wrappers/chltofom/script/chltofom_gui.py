
from qtgui.CCP4TaskWidget import CTaskWidget
from baselayer import QtCore

#-------------------------------------------------------------------
class CTaskchltofom(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'chltofom'
  TASKMODULE = 'expt_data_utility'                               # Where this plugin will appear on the gui
  TASKTITLE = 'Convert HL phase to and from Phi/Fom'             # A short title for gui menu
  DESCRIPTION = '''Interconvert phases between Hendrickson Lattman and Phi/FOM representation'''

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('chltofom')

    self.openFolder(folderFunction='inputData',followFrom=False)

    self.createLine ( [ 'widget','HKLIN' ] )
    
