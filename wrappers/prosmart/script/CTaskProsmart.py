"""
     tasks/prosmart/CTaskProsmart.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class CTaskProsmart(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'prosmart'
  TASKVERSION = 0.1
  TASKMODULE='wrappers'
  TASKTITLE='ProSMART - Restraint generation and structural comparison'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def drawContents(self):

      
    self.setProgramHelpFile('prosmart')
                        
    self.openFolder(folderFunction='inputData')
    self.createLine( [  'widget', '-guiLabel','Target','XYZIN' ] )
    self.createLine( [ 'widget' ,'-guiLabel','Reference', '-browseDb', True, 'REFERENCE_MODEL'] )
    self.createLine( [ 'label', 'External Keywords File: ', 'tip', 'file containing additional prosmart keywords', 'widget', 'EXT_FILE' ] )
    self.createLine( ['widget' , '-guiLabel', 'Fragment library', 'FRAGLIB' ] )
      
    self.openFolder(title='Simple Options')

    self.autoGenerate(container=self.container.controlParameters,selection={'includeParameters' : ['PROGRAM_MODE','LIB_MODE']})
    self.openFolder(title='Advanced options')
  
    self.openTabFrame(title='Alignment')

    self.autoGenerate(container=self.container.controlParameters,selection={'includeParameters' : ['ALIGN_MODE','FRAGLEN','ALIGN_THRESHOLD','HELIX_CUTOFF','HELIX_PENALTY','REWARD_SEQ','ALIGN_REFINE']})
          
    self.openTabFrame(title='Rigid Substructure ID')

    self.autoGenerate(container=self.container.controlParameters,selection={'includeParameters' : [ 'CLUSTER_*' ] } )
    
    self.openTabFrame(title='Restraint Generation')
    self.autoGenerate(container=self.container.controlParameters,selection={'includeParameters' : [ 'RESTRAIN_*' ] } )
    
    self.openTabFrame(title='H-Bond Restraints')
    self.autoGenerate(container=self.container.controlParameters,selection={'includeParameters' : [  'H_*' ] } )
      
    self.openTabFrame(title='Scoring and Output')

    self.autoGenerate(container=self.container.controlParameters,selection={'includeParameters' : ['SUPERPOSE_THRESHOLD','INCLUDE_MAIN','PERFORM_FLIPS','OUTPUT_DM','DISPLAY_AS_DEGREES' ,'OUTPUT_PDB_CHAIN_RESTRAINTS', 'MERGE_CHAINS','RENAME_CHAIN','IS_NMR_MD_ENSEMBLE'] } )

    self.closeTabFrame()

