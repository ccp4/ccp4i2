from __future__ import print_function
#from PyQt4 import QtGui,QtCore
from PySide2 import QtGui,QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

def whatNext(jobId=None):
  return [ 'modelcraft' ]


class Cmrbump_basic(CCP4TaskWidget.CTaskWidget):

  TASKNAME = 'mrbump_basic'
  TASKVERSION = 0.1
  TASKMODULE='molecular_replacement'
  TASKTITLE='Automated structure solution - MrBUMP'
  DESCRIPTION='Run a quick MrBUMP job with streamlined settings'
  RANK=1
  WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']

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

    self.setProgramHelpFile('mrbump_basic')

#-  --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='inputData',title='Input Data')

    self.createLine( [ 'subtitle', 'Target Sequence', 'Sequence of your target structure' ] )

    self.createLine( [ 'widget', 'ASUIN' ] )
    self.createLine( [ 'label', 'The number of monomers to search for', 'widget', 'NMON' ] )

    self.createLine( [ 'subtitle', 'Experimental Data', 'Observed intensities or amplitudes, and Free-R flags' ] )
    self.createLine( [ 'widget', 'F_SIGF' ] )
    self.createLine( [ 'widget', 'FREERFLAG' ] )

#-   --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='modelParameters',title='Search Models')

    #self.createLine( [ 'advice', 'Non-redundancy level for homologue search:' ] )
    #self.createLine( [ 'widget', 'REDUNDANCYLEVEL' ] )
    #self.createLine( [ 'advice', 'Maximum no. of search models to create:' ] )
    #self.createLine( [ 'widget', 'MRMAX' ] )

    # Set up search options
    self.createLine( [ 'subtitle', 'Model databases', 'Databases to search for possible search models' ] )
    self.createLine( [ 'widget', 'SEARCH_PDB', 'label', 'Search PDB for for possible MR search models' ] )
    self.createLine( [ 'advice', indent+'Non-redundancy level for homologue search:' ], toggle=['SEARCH_PDB', 'open', [ True ] ] )
    self.createLine( [ 'label', indent, 'widget', 'REDUNDANCYLEVEL' ], toggle=['SEARCH_PDB', 'open', [ True ] ]  )

    self.createLine( [ 'widget', 'SEARCH_AFDB', 'label', 'Search EBI-AFDB for possible MR search models' ] )
    self.createLine( [ 'advice', indent+'EBI-AFDB pLDDT residue score cut-off:' ], toggle=['SEARCH_AFDB', 'open', [ True ] ]   )
    self.createLine( [ 'label', indent, 'widget', 'AFDBLEVEL' ], toggle=['SEARCH_AFDB', 'open', [ True ] ]   )

    self.createLine( [ 'advice', 'Maximum no. of search models to create:' ] )
    self.createLine( [ 'widget', 'MRMAX' ] )

    self.createLine( [ 'subtitle', 'Optional Settings', 'HHpred results and path to local PDB mirror' ] )
    self.createLine( [ 'widget', '-browseDb', True, 'HHPREDIN', 'tip', 'HHPred results' ] )# ,toggle=['PASTEORREAD','open',['HHPREDIN']])
    self.createLine( [ 'widget', '-browseDb', True, 'PDBLOCAL', 'tip', 'Local PDB mirror' ] )# ,toggle=['PASTEORREAD','open',['HHPREDIN']])

    #self.createLine( [ 'subtitle', 'Local PDB files', 'Locally stored pdb files to be used as MR search models' ] )
    #self.createLine( [ 'widget', 'XYZIN_LIST' ] )

    #Special widget for the input ensemble(s)
    #self.createLine(['subtitle','Search model(s)'])
    #self.openSubFrame(frame=True)
    #self.createLine(['widget', '-title', 'Click "Show list" if more than one copy or more than one search model', 'ENSEMBLES'])
    #self.closeSubFrame()

    #self.createLine( [ 'label','Include specified chains','widget', 'INCLUDE' ] )
    #self.connect ( self.container.controlParameters.INCLUDE, QtCore.SIGNAL('dataChanged'), self.updateRequirements )

    #self.openSubFrame(frame=[True], toggle = ['INCLUDE', 'open', [ True ] ] )
    #self.createLine( ['subtitle','User specified search models'] )
    #self.createLine(['widget','-title','Specify chains to include in MR (e.g. 1smw_A)','SELECTEDCHAINS'])
    #self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Local coordinate files to be used as search models', 'Locally stored coordinate files for inclusion as search models in MR' ] )
    self.createLine( [ 'label','Include local files','widget', 'LOCAL' ] )
    self.openSubFrame(frame=[True], toggle = ['LOCAL', 'open', [ True ] ] )
    #self.createLine( ['subtitle','Local files'] )
    self.createLine ( [ 'tip','Input model from which subset will be selected', 'widget','XYZIN_LIST' ] )
    from qtgui.CCP4ModelWidgets import CPdbDataFileView
    for pdbDataFileView in self.findChildren(CPdbDataFileView):
        pdbDataFileView.showAtomSelection()
    self.createLine( [ 'label','Only use locally provided search models (no sequence search)','widget', 'LOCALONLY' ] )
    self.closeSubFrame()

    self.closeFolder()

    folder = self.openFolder(folderFunction='controlParameters',title='Options')

    self.createLine( [ 'subtitle', 'Molecular Replacement', 'Settings controlling MR step' ] )
    self.createLine( [ 'advice', 'Number of cores for Phaser (maximum=%d)' % MAXPROC ] )
    self.createLine( [ 'widget', 'PJOBS' ] )

    self.createLine( [ 'subtitle', 'Refinement', 'Settings controlling refinement step' ] )
    self.createLine( [ 'advice', 'Number of refinement cycles in Refmac' ] )
    self.createLine( [ 'widget', 'NCYC' ] )

    self.createLine( [ 'subtitle', 'Model Building', 'Settings controlling model building step' ] )
    self.createLine( [ 'advice', 'Run Buccaneer after refinement' ] )
    self.createLine( [ 'widget', 'BUCC' ] )

    self.closeFolder()

#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_PRF', 'advice', 'Search method' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'SAPTF = Spherically Averaged Phased Translation Function' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'PRF = Phased Rotation Function' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'RF(M) = Rotation Function (f-obs from the density outside of the fixed Model)' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'RF(S) = Rotation Function (f-obs from the density inside a Sphere)' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )
#    self.setMenuText( 'PRF', {
#       'n': 'RF(M) + PTF',
#       'y': 'SAPTF + PRF + PTF',
#       's': 'SAPTF + RF(S) + PTF',
#    } )
#    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'PRF' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_PRF', 'open', [ True ] ] )


#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_NPEAKS', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_NPEAKS', 'advice', 'Number of peaks to analyse' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_NPEAKS', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'Number of Rotation Function peaks', 'widget', 'NP'], toggle=[ 'OPEN_NPEAKS', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'Number of Translation Function peaks', 'widget', 'NPT'], toggle=[ 'OPEN_NPEAKS', 'open', [ True ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_NPEAKS', 'open', [ True ] ] )


#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SCORE', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_SCORE', 'advice', 'Scoring putative solutions' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SCORE', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'CC = Correlation Coefficient' ], toggle=[ 'OPEN_SCORE', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'PF = Packing Function' ], toggle=[ 'OPEN_SCORE', 'open', [ True ] ] )
#    self.setMenuText( 'SCORE', {
#       'y': 'Use CC times PF as a score and stop translation search if contrast is > 3.0',
#       'n': 'Use CC times PF and do not stop',
#       'c': 'Use CC and do not stop',
#    } )
#    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'SCORE' ], toggle=[ 'OPEN_SCORE', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'label', 'Expected number of copies (for contrast calculation only)', 'widget', 'NMON_EXP'], toggle=[ 'OPEN_SCORE', 'open', [ True ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SCORE', 'open', [ True ] ] )


#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_ANISO', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_ANISO', 'advice', 'Scaling' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_ANISO', 'open', [ True ] ] )
#    self.setMenuText( 'ANISO', {
#       'y': 'anisotropic',
#       'n': 'isotropic',
#       'k': 'none',
#    } )
#    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'ANISO' ], toggle=[ 'OPEN_ANISO', 'open', [ True ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_ANISO', 'open', [ True ] ] )


#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_HIGH_PATH_VAR', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_HIGH_PATH_VAR', 'advice', 'High pass filter parameter' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_HIGH_PATH_VAR', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'advice', '(B-add, the B-factor applied to input structure amplitudes)' ], toggle=[ 'OPEN_HIGH_PATH_VAR', 'open', [ True ] ] )
#    self.setMenuText( 'HIGH_PATH_VAR', {
#       's': 'From identity between model and sequence (if sequence given)',
#       'i': 'From identity specified manually',
#       'r': 'From high resolution limit',
#       'b': 'Directly as the value of additional B-factor',
#    } )
#    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'HIGH_PATH_VAR' ], toggle=[ 'OPEN_HIGH_PATH_VAR', 'open', [ True ] ] )
#    lab1 = 'Identity between search and target sequences (from 0 to 1)'
#    self.createLine( [ 'label', '         ', 'label', lab1, 'widget', 'SIM' ], toggleFunction=[ self.openSIM, [ 'OPEN_HIGH_PATH_VAR', 'HIGH_PATH_VAR' ] ] )
#    lab1 = 'High resolution limit, Ang'
#    self.createLine( [ 'label', '         ', 'label', lab1, 'widget', 'RESMAX' ], toggleFunction=[ self.openRESMAX, [ 'OPEN_HIGH_PATH_VAR', 'HIGH_PATH_VAR' ] ] )
#    lab1 = 'B-add'
#    self.createLine( [ 'label', '         ', 'label', lab1, 'widget', 'BADD' ], toggleFunction=[ self.openBADD, [ 'OPEN_HIGH_PATH_VAR', 'HIGH_PATH_VAR' ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_HIGH_PATH_VAR', 'open', [ True ] ] )


#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_LOW_PATH_VAR', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_LOW_PATH_VAR', 'advice', 'Low pass filter parameter' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_LOW_PATH_VAR', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'advice', '(B-off, the B-factor of the removed fraction of structure amplitudes)' ], toggle=[ 'OPEN_LOW_PATH_VAR', 'open', [ True ] ] )
#    self.setMenuText( 'LOW_PATH_VAR', {
#       'c': 'From completeness of the search model',
#       'r': 'From low resolution limit',
#       'b': 'Directly as the value of additional B-factor',
#    } )
#    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'LOW_PATH_VAR' ], toggle=[ 'OPEN_LOW_PATH_VAR', 'open', [ True ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_LOW_PATH_VAR', 'open', [ True ] ] )


#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SEQ', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_SEQ', 'advice', 'Sequence modification options' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SEQ', 'open', [ True ] ] )
#    self.createLine( [ 'label', '         ', 'advice', 'Perform alignment and use it to rename residues and trimm side chains' ], toggle=[ 'OPEN_SEQ', 'open', [ True ] ] )
#    self.setMenuText( 'SEQ', {
#       'y': 'always',
#       'd': 'only for sequence identity > 20%',
#       'n': 'never'
#    } )
#    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'SEQ' ], toggle=[ 'OPEN_SEQ', 'open', [ True ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SEQ', 'open', [ True ] ] )


#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SURF', 'open', [ True ] ] )
#    self.createLine( [ 'widget', 'OPEN_SURF', 'advice', 'B-factors modification options' ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SURF', 'open', [ True ] ] )
#    self.setMenuText( 'SURF', {
#       'y': 'Increase B-factor on the molecular surface for all functions',
#       'c': 'Increase B-factor on the surface for Packing function only',
#       'n': 'Do not do anything',
#       '2': 'Set B-factors of all atoms to 20',
#       'a': 'Use poly-alanine model with all B-factors 20',
#    } )
#    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'SURF' ], toggle=[ 'OPEN_SURF', 'open', [ True ] ] )
#-  self.createLine( [ 'advice', '' ], toggle=[ 'OPEN_SURF', 'open', [ True ] ] )


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


