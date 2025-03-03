"""
     tasks/unique/CTaskUnique.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class CTaskUnique(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'unique'
  TASKVERSION = 0.1
  TASKMODULE='test'
  TASKTITLE='Create dummy dataset'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)


  def drawContents(self):

    self.setProgramHelpFile('unique')
    
    self.openFolder(folderFunction='inputData')

    self.createLine(  [  'advice', 'This task will create a dummy list of reflections'] )

    self.createLine ( [  'widget','SPACEGROUPCELL' ] )

    self.createLine ( [  'label', 'Generate reflections to resolution:',
                         'widget','RESOLUTION' ] )
