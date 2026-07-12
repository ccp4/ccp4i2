from qtgui.CCP4TaskWidget import CTaskWidget

class Ccoot_find_waters(CTaskWidget):
  TASKNAME = 'coot_find_waters'
  TASKVERSION = 0.1
  TASKMODULE='model_building'
  TASKTITLE='Find waters - COOT'
  DESCRIPTION='Find and filter waters based on electron density and contacts (non-interactive Coot)'

  def drawContents(self):
    self.setProgramHelpFile('coot_find_waters')
    self.openFolder(folderFunction='inputData')

    self.createLine(["subtitle", "Input data"])
    self.openSubFrame(frame=True)
    self.createLine(['tip', 'Input structure', 'widget', 'XYZIN'])
    self.createLine(['tip', 'Input 2mFo-DFc coefficients', 'widget', 'FPHIIN'])
    self.closeSubFrame()

    self.createLine(["subtitle", "Options"])
    self.openSubFrame(frame=True)
    self.createLine(['label', 'Threshold: ', 'tip', 'Threshold for water addition (in units of RMS)', 'stretch','widget', 'THRESHOLD'])
    self.createLine(['label', 'Minimum distance: ', 'tip', 'Minimum distance from existing atoms', 'stretch','widget', 'MINDIST'])
    self.createLine(['label', 'Maximum distance: ', 'tip', 'Maximum distance from existing atoms', 'stretch','widget', 'MAXDIST'])
    self.closeSubFrame()
