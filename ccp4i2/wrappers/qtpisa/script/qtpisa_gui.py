"""
     qtpisa task widget
"""

from ....qtgui.CCP4TaskWidget import CTaskWidget


class Cqtpisa(CTaskWidget):

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
    #print 'Ccoot_rebuild.isValid'
    if self.getWidget('followFrom') is None: return
    followJobId = self.getWidget('followFrom').currentJobId()
    #print 'Ccoot_rebuild.isValid followFrom',followJobId
    return CTaskWidget.isValid(self)
