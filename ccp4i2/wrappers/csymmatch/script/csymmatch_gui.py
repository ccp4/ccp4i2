

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class csymmatch_gui(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'csymmatch'
  TASKVERSION = 0.0
  TASKMODULE=[ 'molecular_replacement', 'model_data_utility' ]
  TASKTITLE='Match model to reference structure'
  SHORTTASKTITLE='MATCH MODEL'
  DESCRIPTION='Match symmetry and origin of output model to reference structure (Csymmatch)'
  WHATNEXT = []
  PROGRAMHELP = 'csymmatch'
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)


  def drawContents(self):

    self.setProgramHelpFile('csymmatch')

    
    # Remove the 'followFrom' widget cos there can be no preceeding jobs - (this may be wrong thing to do)                    
    folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)

    self.createLine( [ 'subtitle', 'Atomic model to be moved' ] )
    self.createLine( ['widget', 'XYZIN_QUERY'] )

    self.createLine( [ 'subtitle', 'Fixed (reference) atomic model' ] )
    self.createLine( ['widget', 'XYZIN_TARGET'] )

    self.createLine( [ 'label','Try all possible origins and hands','widget', 'ORIGIN_HAND' ] )
    self.createLine( [ 'label','Radius to use in stiching floating fragments to chains','widget', 'CONNECTIVITY_RADIUS' ] )
