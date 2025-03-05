from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Modules

class Ccoot1(CTaskWidget):

  TASKNAME = 'coot1'
  TASKVERSION = 0.1
  TASKMODULE='model_building'
  TASKTITLE='Coot 1'
  SHORTTASKTITLE='COOT1'
  DESCRIPTION='Interactive building (Coot)'

  WHATNEXT = ['prosmart_refmac','coot_rebuild', 'modelcraft' ]

  def drawContents(self):
    
    self.setProgramHelpFile('coot_rebuild')
                        
    self.openFolder(folderFunction='inputData')
    self.createLine(['subtitle','Key bindings'])
    self.openSubFrame(frame=True)
    self.createLine( [ 'label', 'Use Bernhard and Paul Key bindings','stretch','widget', 'USEKEYBINDINGS' ] )
    self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Coordinates' ] )
    self.createLine( [ 'widget', 'XYZIN_LIST' ] )
    self.createLine( [ 'subtitle', 'Electron density maps' ] )
    self.createLine( [ 'widget', 'FPHIIN_LIST' ] )
    self.createLine( [ 'subtitle', 'Difference density maps' ] )
    self.createLine( [ 'widget', 'DELFPHIIN_LIST' ] )
    self.createLine( [ 'subtitle', 'Anomalous difference density maps' ] )
    self.createLine( [ 'widget', 'DELFPHIINANOM_LIST' ] )
    self.createLine( [ 'subtitle', 'Additional data' ] )
    self.createLine( [ 'widget', 'DICT' ] )
    self.createLine( [ 'widget', 'COOTSCRIPTFILE' ] )
