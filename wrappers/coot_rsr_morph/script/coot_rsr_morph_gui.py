from qtgui.CCP4TaskWidget import CTaskWidget

class Ccoot_rsr_morph(CTaskWidget):
  TASKNAME = 'coot_rsr_morph'
  TASKVERSION = 202110070909
  TASKMODULE='refinement'
  TASKTITLE='Real space refinement morphing with coot'
  DESCRIPTION='Real space refinement morphing with coot'

  def drawContents(self):
    self.setProgramHelpFile('coot_rsr_morph')
    self.openFolder(folderFunction='inputData')
    self.createLine( [ 'subtitle', 'Main inputs' ])
    self.openSubFrame(frame=[True])
    self.createLine( [ 'tip', 'input structure', 'widget', 'XYZIN' ] )
    self.createLine( [ 'tip', 'input 2mFo-DFc coefficient', 'widget', 'FPHIIN' ] )
    self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Options' ])
    self.openSubFrame(frame=[True])
    self.createLine( [ 'label', 'Local radius', 'widget', 'LOCAL_RADIUS' ] )
    self.createLine( [ 'label', 'GM alpha', 'widget', 'GM_ALPHA' ] )
    self.createLine( [ 'label', 'Blur B-factor', 'widget', 'BLUR_B_FACTOR' ] )
    self.closeSubFrame()
