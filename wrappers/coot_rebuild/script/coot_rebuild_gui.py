"""
     coot_rebuild task widget
"""

from baselayer import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Modules

#-------------------------------------------------------------------
class Ccoot_rebuild(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'coot_rebuild'
  TASKVERSION = 0.1
  TASKMODULE='model_building'
  TASKTITLE='Model building - COOT 0.9'
  SHORTTASKTITLE='COOT'
  DESCRIPTION='Interactive building (Coot 0.9)'

  WHATNEXT = ['prosmart_refmac','coot_rebuild', 'modelcraft', 'servalcat_pipe' ]

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

  
  def isValid(self):
    import os
    #print 'Ccoot_rebuild.isValid'
    if self.getWidget('followFrom') is None: return
    followJobId = self.getWidget('followFrom').currentJobId()
    #print 'Ccoot_rebuild.isValid followFrom',followJobId
    if followJobId  is not None:
      stateFile = os.path.join( CCP4Modules.PROJECTSMANAGER().db().jobDirectory(jobId=followJobId),'COOT_FILE_DROP','0-coot.state.scm')
      if os.path.exists(stateFile):
        self.container.inputData.COOTSTATEFILE.setFullPath(stateFile)
    return CTaskWidget.isValid(self)
  
