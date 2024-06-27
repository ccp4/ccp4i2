from __future__ import print_function
#from PyQt4 import QtGui,QtCore
from PySide2 import QtGui,QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

def whatNext(jobId=None):
  return [ 'buccaneer_build_refine_mr' ]


class Cslicendice(CCP4TaskWidget.CTaskWidget):

  TASKNAME = 'slicendice'
  TASKVERSION = 0.1
  TASKMODULE=['alpha_fold']
  TASKTITLE='SliceNDice - Auto model processing and MR'
  DESCRIPTION='Automated processing of predicted or deposited search models and Molecular Replacement'
  RANK=1
  WHATNEXT = ['prosmart_refmac','buccaneer_build_refine_mr','coot_rebuild']

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def updateRequirements ( self ) :

    if self.container.controlParameters.INCLUDE :
      self.container.inputData.SELECTEDCHAINS.setQualifier ( 'allowUndefined', False )
      self.container.inputData.SELECTEDCHAINS.setQualifier ( 'toolTip', 'A chain is required' )
    else :
      self.container.inputData.SELECTEDCHAINS.setQualifier ( 'allowUndefined', True )
    self.getWidget ( 'SELECTEDCHAINS' ).validate ( )
    
    return

  def drawContents(self):

    import multiprocessing
    MAXPROC=multiprocessing.cpu_count()  
 
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

    self.setProgramHelpFile('slicendice')

#-  --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='inputData',title='Input Data')

    self.createLine( ['advice', '<b>SliceNDice is an automated pipeline for processing models for use in Molecular Replacement (MR). It can correct B-factor \
           <br> columns in predicted models as well "slicing" search models into rigid regions that may achieve better placement in MR. \
           <br> Following this, it will attempt MR using all models created. </b>'] )

    self.createLine( [ 'subtitle', 'Experimental Data', 'Observed intensities or amplitudes, and Free-R flags' ] )
    self.createLine( [ 'widget', 'F_SIGF' ] )

    self.createLine( [ 'widget', 'FREERFLAG' ] )

    self.createLine( [ 'subtitle', 'Target Sequence', 'Sequence of your target structure' ] )

    self.createLine( [ 'widget', 'ASUIN' ] )
    self.createLine( [ 'label', 'The number of monomers to search for', 'widget', 'NO_MOLS' ] )

    self.createLine( [ 'subtitle', 'Search Model', 'Molecular Replacement search model (e.g. Alphafold or RosettaFold prediction, deposited pdb structure)' ] )

    self.createLine( [ 'widget', 'XYZIN' ] )

    self.createLine(['advice', 'Select B-factor treatment option - it is important this is set correctly' ])
    self.setMenuText('BFACTOR_TREATMENT', {'plddt': 'Convert pLDDT scores to B-factors (typical of Alphafold predictions)',
                                           'rms': 'Convert rmsd estimates to B-factors (RosettaFold)',
                                           'fractional_plddt': 'Convert fractional pLDDT scores to B-factors (ESMFold predicitons)',
                                           'predicted_bfactor': 'Predicted B-factors (no conversion required)',
                                           'bfactor': 'B-factors already present (e.g. deposited PDB structure)' })
    self.createLine(['label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'BFACTOR_TREATMENT']) 

    self.createLine( [ 'advice', 'pLDDT residue score threshold:' ] , toggle=['BFACTOR_TREATMENT', 'open', [ 'plddt' ] ]   )
    self.createLine( [ 'label', indent, 'widget', 'PLDDT_THRESHOLD' ] , toggle=['BFACTOR_TREATMENT', 'open', [ 'plddt' ] ]   )
    self.createLine( [ 'advice', 'rms residue score threshold:' ] , toggle=['BFACTOR_TREATMENT', 'open', [ 'rms' ] ]   )
    self.createLine( [ 'label', indent, 'widget', 'RMS_THRESHOLD' ] , toggle=['BFACTOR_TREATMENT', 'open', [ 'rms' ] ]   )

#-   --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='controlParameters',title='Options')

    self.createLine( [ 'subtitle', 'Model slicing', 'Settings controlling slicing of search models' ] )

    self.createLine( ['advice', 'These parameters control model splitting. MIN and MAX splits are the minimum and maximum number of splits \
                                 <br>to make to the search models. For example, a MIN split of 1 and a MAX split of 3 will result in 3 MR jobs. \
                                 <br>(1) using an unsplit model, (2) a model split in two and (3) a model split three ways.'] )

    self.createLine( [ 'advice', 'Minimum split level' ] )
    self.createLine( [ 'widget', 'MIN_SPLITS' ] )

    self.createLine( [ 'advice', 'Maximum split level' ] )
    self.createLine( [ 'widget', 'MAX_SPLITS' ] )

    self.createLine( [ 'subtitle', 'Molecular Replacement', 'Settings controlling MR step' ] )
    self.createLine( [ 'advice', 'Number of cores for Phaser (maximum=%d)' % MAXPROC ] )
    self.createLine( [ 'widget', 'NPROC' ] )

    self.createLine( [ 'subtitle', 'Refinement', 'Settings controlling refinement step' ] )
    self.createLine( [ 'advice', 'Number of refinement cycles in Refmac' ] )
    self.createLine( [ 'widget', 'NCYC' ] )

    self.closeFolder()

#-   --------------------          --------------------          --------------------

    print('CTaskMolrep stackedWidgets')
    for w in self.findChildren( CCP4TaskWidget.CStackedWidget ) :
      print('   ', w, w.controlVar)

    self.updateViewFromModel()

#-   --------------------          --------------------          --------------------


  def openSIM( self ) :

    gui = self.container.guiParameters
    if gui.OPEN_HIGH_PATH_VAR and str( gui.HIGH_PATH_VAR ) == 'i' :
      return True

    else:
      return False


  def openRESMAX( self ) :

    gui = self.container.guiParameters
    if gui.OPEN_HIGH_PATH_VAR and str( gui.HIGH_PATH_VAR ) == 'r' :
      return True

    else:
      return False


  def openBADD( self ) :

    gui = self.container.guiParameters
    if gui.OPEN_HIGH_PATH_VAR and str( gui.HIGH_PATH_VAR ) == 'b' :
      return True

    else:
      return False


