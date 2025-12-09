from ccp4i2.baselayer import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets
from ccp4i2.core import CCP4Modules

class CTaskPDB_REDO(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'pdb_redo_api'
  TASKVERSION = 0.0
  TASKMODULE='refinement'
  TASKTITLE='PDB-REDO Web services'
  SHORTTASKTITLE='PDB-REDO'
  WHATNEXT = []
  PROGRAMHELP = 'pdb_redo'
  DESCRIPTION='Refine structures using PDB-REDO Web service'
  
  def drawContents(self):

    self.setProgramHelpFile('pdb_redo')
    
    folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)

    if not CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_SECRET or not CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID or not CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID:
        self.createLine( [ 'subtitle', '<span style="color:red;">You do not have valid PDB-REDO token/secret set. Please set these using:<br/> <em>Utilities -> System administration tools -> Configure login tokens for PDB-REDO.<br/><br/>You can obtain a token/secret pair from https://services.pdb-redo.eu/</em></span>' ] )

    self.createLine( [ 'subtitle', 'Main Input data' ] )
    self.openSubFrame(frame=[True])
    self.createLine( ['tip', 'This is the structure which will be refined.', 'widget', 'XYZIN'] )
    self.createLine( [ 'widget', 'F_SIGF' ] )
    self.createLine( [ 'widget', 'FREERFLAG' ] )
    self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Optional input data' ] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'SEQIN' ] )
    self.createLine( [ 'widget', 'DICT' ] )
    self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Paired refinement' ] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'PAIRED', 'label', 'Perform paired refinement' ] )
    self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Rebuilding Options' ] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'NOLOOPS', 'label', 'Do not try to complete loops' ] )
    self.createLine( [ 'widget', 'NOPEPFLIP', 'label', 'Do not perform peptide flips' ] )
    self.createLine( [ 'widget', 'NOSCBUILD', 'label', 'Do not rebuild or complete side-chain' ] )
    self.createLine( [ 'widget', 'NOCENTRIFUGE', 'label', 'Do not delete poor waters' ] )
    self.createLine( [ 'widget', 'NOSUGARBUILD', 'label', 'Do not (re)build carbohydrates' ] )
    self.createLine( [ 'widget', 'NOREBUILD', 'label', 'Skip model rebuilding' ] )
    self.closeSubFrame()

    folder = self.openFolder(folderFunction='controlParameters',title='Advanced Options')
    self.createLine( [ 'subtitle', 'Advanced Options' ] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'NEWMODEL', 'label', 'Update the model even if the initial refinement is poor' ] )
    self.createLine( [ 'widget', 'ISOTROPIC', 'label', 'Force isotropic B-factors' ] )
    self.createLine( [ 'widget', 'ANISOTROPIC', 'label', 'Force anisotropic B-factors (within limits)' ] )
    self.createLine( [ 'widget', 'NOTLS', 'label', 'Do not perform TLS refinement' ] )
    self.createLine( [ 'widget', 'TIGHTER', 'label', 'Use tighter restraints in refinement' ] )
    self.createLine( [ 'widget', 'LOOSER', 'label', 'Use looser restraints in refinement' ] )
    self.createLine( [ 'advice', 'Note: The tighter and looser values are subtracted so the net effect of having both set to 2 will be none' ] )
    self.closeSubFrame()
 
  def isValid(self):
    if not CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_SECRET or not CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID or not CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID:

        invalidElements = super(CTaskPDB_REDO, self).isValid()
        invalidElements.append(CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_SECRET)

        return invalidElements
        
    return CCP4TaskWidget.CTaskWidget.isValid(self)
 
