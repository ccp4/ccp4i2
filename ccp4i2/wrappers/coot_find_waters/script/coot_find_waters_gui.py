from ......qtgui.CCP4TaskWidget import CTaskWidget


#-------------------------------------------------------------------
class Ccoot_find_waters(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'coot_find_waters'
  TASKVERSION = 0.1
  TASKMODULE='model_building'
  TASKTITLE='Find waters - COOT'
  DESCRIPTION='Find and filter waters based on electron density and contacts (non-interactive Coot)'

  def drawContents(self):

    self.setProgramHelpFile('coot_find_waters')

    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'tip', 'input structure', 'widget', 'XYZIN' ] )
    self.createLine( [ 'tip', 'input 2mFo-DFc coefficient', 'widget', 'FPHIIN' ] )
 
    self.openFolder(folderFunction='controlParameters',title='Options')
    self.createLine( [ 'label', 'Threshold: ', 'tip', 'Threshold for water addition(in units of RMS)', 'stretch','widget', 'THRESHOLD' ] )
    self.createLine( [ 'label', 'Minimum distance: ', 'tip', 'Minimum distance from existing atom', 'stretch','widget', 'MINDIST' ] )
    self.createLine( [ 'label', 'Maximum distance: ', 'tip', 'Minimum distance from existing atom', 'stretch','widget', 'MAXDIST' ] )
