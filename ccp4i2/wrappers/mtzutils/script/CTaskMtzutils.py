"""
     tasks/mtzutils/CTaskMtzutils.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class CTaskMtzutils(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'mtzutils'
  TASKVERSION = 0.0
  TASKMODULE='test'
  TASKTITLE='Add or delete MTZ columns'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)


  def drawContents(self):

    self.setProgramHelpFile('mtzutils')

    self.openFolder(folderFunction='protocol')

    self.createTitleLine()

    self.createLine( [ 'widget', 'MODE',
                       'label', 'from an MTZ file'
                       ] )
    
    self.openFolder(folderFunction='inputData')

    self.createLine ( [  'help', 'TESTING',
                         'widget','HKLIN1' ] )

    self.createLine(  [ 'tip','Columns to be included/excluded',
                        'widget','FSIGF'] )

    self.createLine ( [  'help', 'TESTING',
                         'widget','HKLIN2' ] )

    self.openFolder( folderFunction='outputData')

    self.createLine ( [  'help', 'TESTING',
                         'widget','HKLOUT' ] )

