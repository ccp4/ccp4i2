from ....qtgui.CCP4TaskWidget import CTaskWidget


class Cqtpisa(CTaskWidget):
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
    if self.getWidget('followFrom') is None:
      return
    self.getWidget('followFrom').currentJobId()
    return CTaskWidget.isValid(self)
