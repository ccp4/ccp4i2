

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets
import qtgui

class CTaskGesamt(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'gesamt'
  TASKVERSION = 0.0
  TASKMODULE='model_data_utility'
  TASKTITLE='Structural alignment - Gesamt'
  SHORTTASKTITLE='GESAMT'
  WHATNEXT = []
  PROGRAMHELP = 'gesamt'
  DESCRIPTION='Superpose one protein structure on another'
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('gesamt')

    
    # Remove the 'followFrom' widget cos there can be no preceeding jobs - (this may be wrong thing to do)                    
    folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)

    self.createLine( [ 'subtitle', 'Pairwise (2 structures) or multiple alignment' ] )
    self.createLine( ['tip', 'Pair wise or multiple alignment', 'widget', 'PAIRMULTI'] )
    self.openSubFrame(frame=[True],toggle=["PAIRMULTI", "open", ["MULTIPLE"]])
    self.createLine( [ 'subtitle', 'Structures to superpose' ] )
    listWidget = self.createLine( ['tip', 'The structures willbe transformed to a consensus position.', 'widget', 'XYZIN_LIST'] )

    #List widgets are odd. I cannot get them to behave properly - just making white.
    childer = listWidget.findChildren(QtWidgets.QWidget)
    for child in childer:
        if type(child) is qtgui.CCP4Widgets.CListViewListWidget:
            child.setStyleSheet("QFrame { background-color:white}")

    self.closeSubFrame()
    self.openSubFrame(frame=[True],toggle=["PAIRMULTI", "open", ["PAIR"]])
    self.createLine( [ 'subtitle', 'Structure to move' ] )
    self.createLine( ['tip', 'This is the query structure to which the transformation matrix will be applied.', 'widget', 'XYZIN_QUERY'] )
    atWidget = self.getWidget('XYZIN_QUERY')
    if hasattr(atWidget, 'showAtomSelection'): self.getWidget('XYZIN_QUERY').showAtomSelection()

    self.createLine( [ 'subtitle', 'Fixed structure' ] )
    self.createLine( ['tip', 'This is the target structure onto which the query structure will be moved.', 'widget', 'XYZIN_TARGET'] )
    atWidget = self.getWidget('XYZIN_TARGET')
    if hasattr(atWidget, 'showAtomSelection'): self.getWidget('XYZIN_TARGET').showAtomSelection()
    self.closeSubFrame()

    #self.createLine( ['tip', 'Save the moved model, or just report the quality of the match?', 'widget', 'OUTPUT_COORDS', 'label', 'Output superposed coordinates'])

    self.openSubFrame(frame=[True])
    self.createLine( ['tip', 'Normal mode is sufficient for most purposes. High quality mode will take about 10 times longer', 'widget', 'MODE' ] )
    self.closeSubFrame()

  def isValid(self):
      if str(self.container.controlParameters.PAIRMULTI) == "MULTIPLE":
          invalidElements = [self.container.inputData.XYZIN_LIST]
          if self.container.inputData.XYZIN_LIST.isSet() and len(self.container.inputData.XYZIN_LIST)>1:
              allSet = True
              for f in self.container.inputData.XYZIN_LIST:
                  allSet = allSet and f.isSet()
              print("alllSet?",allSet)
              if allSet:
                  invalidElements = []
      else:
          invalidElements = []
          if not self.container.inputData.XYZIN_TARGET.isSet():
              invalidElements.append(self.container.inputData.XYZIN_TARGET)
          if not self.container.inputData.XYZIN_QUERY.isSet():
              invalidElements.append(self.container.inputData.XYZIN_QUERY)
      return invalidElements
