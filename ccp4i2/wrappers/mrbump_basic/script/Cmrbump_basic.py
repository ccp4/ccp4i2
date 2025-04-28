import multiprocessing

from ....qtgui import CCP4TaskWidget


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

    self.createLine( [ 'subtitle', 'Local coordinate files to be used as search models', 'Locally stored coordinate files for inclusion as search models in MR' ] )
    self.createLine( [ 'label','Include local files','widget', 'LOCAL' ] )
    self.openSubFrame(frame=[True], toggle = ['LOCAL', 'open', [ True ] ] )
    #self.createLine( ['subtitle','Local files'] )
    self.createLine ( [ 'tip','Input model from which subset will be selected', 'widget','XYZIN_LIST' ] )
    from ....qtgui.CCP4ModelWidgets import CPdbDataFileView
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

    print('CTaskMolrep stackedWidgets')
    for w in self.findChildren( CCP4TaskWidget.CStackedWidget ) :
      print('   ', w, w.controlVar)
    self.updateViewFromModel()

  def openSIM(self):
    gui = self.container.guiParameters
    if gui.OPEN_HIGH_PATH_VAR and str( gui.HIGH_PATH_VAR ) == 'i' :
      return True
    else:
      return False


  def openRESMAX(self):
    gui = self.container.guiParameters
    if gui.OPEN_HIGH_PATH_VAR and str( gui.HIGH_PATH_VAR ) == 'r' :
      return True
    else:
      return False


  def openBADD(self):
    gui = self.container.guiParameters
    if gui.OPEN_HIGH_PATH_VAR and str( gui.HIGH_PATH_VAR ) == 'b' :
      return True
    else:
      return False
