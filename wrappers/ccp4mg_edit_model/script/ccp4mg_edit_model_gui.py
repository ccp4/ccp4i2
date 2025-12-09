"""
     ccp4mg_edit_model task widget
"""

from ccp4i2.baselayer import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.core import CCP4Modules

#-------------------------------------------------------------------
class Cccp4mg_edit_model(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'ccp4mg_edit_model'
  TASKVERSION = 0.1
  TASKMODULE=['alpha_fold','bioinformatics']
  TASKTITLE='Interactive model preparation - CCP4mg and MrBUMP'
  SHORTTASKTITLE='CCP4mg MrBUMP'
  DESCRIPTION='Identify MR models with MrBUMP, display and select with CCP4mg'
  WHATNEXT = ['phaser_simple', 'phaser_pipeline', 'molrep_pipe']

  def drawContents(self):
    
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

    self.setProgramHelpFile('ccp4mg_edit_model')
                        
    self.openFolder(folderFunction='inputData')

    self.createLine( [ 'advice', 'Sequences from AU content:' ] )
    self.createLine( [ 'widget', 'ASUIN' ] )

    # Set up search options
    self.createLine( [ 'subtitle', 'Model databases', 'Databases to search for possible search models' ] )
    self.createLine( [ 'widget', 'SEARCH_PDB', 'label', 'Search PDB for for possible MR search models' ] )
    self.createLine( [ 'advice', indent+'Non-redundancy level for homologue search:' ], toggle=['SEARCH_PDB', 'open', [ True ] ] )
    self.createLine( [ 'label', indent, 'widget', 'REDUNDANCYLEVEL' ], toggle=['SEARCH_PDB', 'open', [ True ] ]  )

    self.createLine( [ 'widget', 'SEARCH_AFDB', 'label', 'Search EBI-AFDB for possible MR search models' ] )
    self.createLine( [ 'advice', indent+'EBI-AFDB pLDDT residue score cut-off:' ], toggle=['SEARCH_AFDB', 'open', [ True ] ]   )
    self.createLine( [ 'label', indent, 'widget', 'AFDBLEVEL' ], toggle=['SEARCH_AFDB', 'open', [ True ] ]   )

    #self.createLine( [ 'advice', 'Cutoff threshold for phmmer results:' ] )
    #self.createLine( [ 'widget', 'PHMMERCUTOFF' ] )
    self.createLine( [ 'advice', 'Maximum no. of search models to create:' ] )
    self.createLine( [ 'widget', 'MRMAX' ] )

    self.createLine( [ 'subtitle', 'Optional Settings', 'HHpred results and path to local PDB mirror' ] )
    self.createLine( [ 'widget', '-browseDb', True, 'HHPREDIN', 'tip', 'HHPred results' ] )# ,toggle=['PASTEORREAD','open',['HHPREDIN']])
    self.createLine( [ 'widget', '-browseDb', True, 'PDBLOCAL', 'tip', 'Local PDB mirror' ] )# ,toggle=['PASTEORREAD','open',['HHPREDIN']])

  def ToggleSPDB(self):
    return str(self.container.inputData.REDUNDANCYLEVEL) == '100'
 
  def isValid(self):
    import os
    #print 'Ccoot_rebuild.isValid'
    if self.getWidget('followFrom') is None: return
    followJobId = self.getWidget('followFrom').currentJobId()
    #print 'Ccoot_rebuild.isValid followFrom',followJobId
    return CTaskWidget.isValid(self)
  
