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

  def updateXYZ2(self):
        if self.container.controlParameters.DIFFMODE.isSet():
            if str(self.container.controlParameters.DIFFMODE) != "COMPARE":
                self.container.inputData.XYZIN2.unSet()

  def drawContents(self):

    self.setProgramHelpFile('areaimol')

    folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)

    self.createLine(['subtitle','Input'])
    self.openSubFrame( frame=[True])

    self.createLine( [ 'subtitle', 'Input Structure' ] )
    self.createLine( ['tip', 'This is the structure which the area will be calculated.', 'widget', 'XYZIN'] )
    self.createLine( ['tip', 'This keyword controls the program function, the data required, and how it is processed and analysed', 'widget', 'DIFFMODE'] )
    self.createLine( ['label', 'The symmetry of the molecule, e.g. P212121', 'widget', 'SYMMETRY'] , toggle=[ 'DIFFMODE', 'open' , ['IMOL'] ] )
    self.createLine( [ 'subtitle', 'Second structure to compare with' ] , toggle=[ 'DIFFMODE', 'open' , ['COMPARE'] ])
    self.createLine( ['tip', 'A second set of input coordinates, and is only used in DIFFMODE COMPARE.', 'widget', 'XYZIN2'] , toggle=[ 'DIFFMODE', 'open' , ['COMPARE'] ] )
    self.createLine( ['label', 'Output mode', 'widget', 'OUTPUT_MODE'], toggle=[ 'DIFFMODE', 'open' , ['OFF'] ]  )
    self.createLine( ['label', 'Output mode', 'widget', 'OUTPUT_MODE_COMPARE'], toggle=[ 'DIFFMODE', 'open' , ['COMPARE','IMOL'] ]  )
    self.closeSubFrame()

    self.connectDataChanged('DIFFMODE', self.updateXYZ2)

    #self.createLine( [ 'subtitle', 'Keywords' ] )
    #self.createLine( ['tip', 'Script'] )
    #self.createLine( ['widget', '-guiMode','multiLine', 'EXTRA_AREAIMOL_KEYWORDS'] )
