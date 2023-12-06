
from qtgui.CCP4TaskWidget import CTaskWidget

class CRefmacWrk(CTaskWidget):

  TASKNAME = 'refmac_wrk'
  TASKVERSION = 0.0
  TASKMODULE = 'test'
  TASKTITLE = 'Refmac Replica'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def showFP( self ) :
    inp = self.container.inputData
    par = self.container.controlParameters
    return str(par.INPUT_PHASE) in ['HL','PHASE'] \
        or str(par.INPUT_PHASE) == 'NO' and str(par.TWINREF_TYPE) in ['NO','AMPLITUDES']

  def showIP( self ) :
    inp = self.container.inputData
    par = self.container.controlParameters
    return str(par.INPUT_PHASE) == 'NO' and str(par.TWINREF_TYPE) == 'INTENSITIES'

  def drawContents(self):

# ----------------------------------------------------------------------

    self.openFolder(
        folderFunction='inputData',
        title='Input Data')

    self.createLine([
        'label','Job title',
        'widget','TITLE'])

    self.createLine([
        'label','Do',
        'widget','REFINE_TYPE',
        'label','using',
        'widget','INPUT_PHASE',
        'label','input'],
        toggle=['REFINE_TYPE','hide',['REVIEW','IDEA']])

    self.createLine([
        'label','Do',
        'widget','REFINE_TYPE'],
        toggle=['REFINE_TYPE','open',['REVIEW','IDEA']])

    self.createLine([
        'widget','IFFIXTLS',
        'label','Input fixed TLS parameters'],
        toggle=['REFINE_TYPE','hide',['REVIEW','IDEA','TLS']])

    self.openSubFrame(
        toggle=['INPUT_PHASE','open',['NO']])

    self.createLine([
        'widget','TWINREF_TYPE',
        'label','twin refinement'])

    self.closeSubFrame()

    self.openSubFrame(
        toggle=['PROSMART_AVAIL','open',[True]])

    self.createLine([
        'widget','IFPROSMART',
        'label','Run Prosmart to generate ',
        'widget','PROSMART_MODE',
        'label','(low resolution refinement)'],
        toggle=['REFINE_TYPE','open',['REST','TLS']])

    self.closeSubFrame()

    self.openSubFrame(
        toggle=['COOT_AVAIL','open',[True]])

    self.createLine([
        'widget','RUN_COOT_FW',
        'label','Run Coot:findwaters to automatically add/remove waters to refined structure'],
        toggle=['REFINE_TYPE','open',['REST','UNRE']])

    self.closeSubFrame()

    self.openSubFrame(
        toggle=['REFINE_TYPE','hide',['REVIEW','IDEA']])

    self.createLine([
        'label','MTZ in',
        'widget','HKLIN'])

    self.createLine([
        'label','F',
        'widget','F_ANO'],
        toggle=['INPUT_PHASE','open',['SAD']])

    self.createLine([
        'label','Wavelength ',
        'widget','WAVELENGTH',
        'label','Refine substructure occupancy ',
        'widget','REF_SUBOCC'],
        toggle=['INPUT_PHASE','open',['SAD']])

    self.createLine([
        'widget', 'ANOMALOUS_ATOMS',
        '-title','Anomalous form factors'],
        toggle=['INPUT_PHASE','open',['SAD']])

    self.createLine([
        'label','IP    ',
        'widget','IOBS'],
        toggleFunction=[self.showIP,['INPUT_PHASE','TWINREF_TYPE']])

    self.createLine([
        'label','FP    ',
        'widget','FOBS'],
        toggleFunction=[self.showFP,['INPUT_PHASE','TWINREF_TYPE']])

    self.createLine([
        'label','PHIB  ',
        'widget','PHOBS'],
        toggle=['INPUT_PHASE','open',['PHASE']])

    self.createLine([
        'label','HLA',
        'widget','ABCD'],
        toggle=['INPUT_PHASE','open',['HL']])

    self.closeSubFrame()

    self.createLine([
        'label','PDB in',
        'widget','XYZIN'])

    self.createLine([
        'label','LIB in',
        'widget','MAKE_LIBRARY'],
        toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])

#  add extra button to launch merge_monomers task

    self.openSubFrame(
        toggle=['REFINE_TYPE','open',['TLS']])

    self.createLine([
        'label','TLS in (optional) ',
        'widget','TLSIN'])

#  add extra button to launch merge_monomers task

    self.closeSubFrame()

    self.openSubFrame(
        toggle=['IFFIXTLS','open',[True]])

    self.createLine([
        'label','TLS in',
        'widget','TLSIN'])

    self.closeSubFrame()

    self.openSubFrame(
        toggle=['IFPROSMART','open',[True]])

    self.createLine([
        'label','Reference PDB in',
        'widget','EXT_XYZIN'],
        toggle=['PROSMART_MODE','hide',['SECSTR']])

    self.createLine([
        'label','Prosmart restraints out',
        'widget','RESTRAINTFILE'])

    self.createLine([
        'label','Prosmart keyword file',
        'widget','PROSMART_KEYFILE'])

    self.closeSubFrame()

    self.openSubFrame()

    self.createLine([
        'label','Refmac keyword file',
        'widget','INCLUDEFILE'])

    self.closeSubFrame()


# ----------------------------------------------------------------------

    self.openFolder(
        folderFunction='controlParameters',
        title='Options')

#   self.drawDataHarvesting()
    self.drawSetupGeometricRestraints()
    self.drawSetupGeometricRestraintsReview()
    self.drawSetupNCSRestraints()
    self.drawRefinementParametersREST()
    self.drawRefinementParametersUNRE()
    self.drawRefinementParametersRIGID()
    self.drawRefinementParametersTLS()
    self.drawExternalRestraints()
    self.drawIdealisationParameters()
    self.drawCootParameters()
    self.drawTLSParameters()
    self.drawMonitoringAndOutputOptions()
    self.drawScaling()
    self.drawRigidDomainsDefinition()
    self.drawGeometricParameters()

    self.openSubFrame(
        frame=True)
    self.closeSubFrame()

# ----------------------------------------------------------------------
# def drawDataHarvesting(self):

#   self.openSubFrame(
#       frame=True,
#       toggle=['REFINE_TYPE','hide',['REVIEW','RIGID','IDEA']])
#   self.closeSubFrame()

#   self.openSubFrame(
#       frame=True,
#       toggle=['REFINE_TYPE','hide',['REVIEW','RIGID','IDEA']])

#   self.createLine([
#       'widget','DATA_HARVESTING',
#       'advice','Data Harvesting'])

#   self.createLine([
#       'label','AAAAAAAAAAAAAAAAAAAAAA 01'],
#       toggle=['DATA_HARVESTING','open',[True]])

#   self.closeSubFrame()

# ----------------------------------------------------------------------
  def drawSetupGeometricRestraints(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['REVIEW','RIGID','UNRE']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['REVIEW','RIGID','UNRE']])

    self.createLine([
        'widget','SETUP_GEOMETRIC_RESTRAINTS',
        'advice','Setup Geometric Restraints'])

    self.createLine([
        'label','Checking against dictionary:',
        'widget','MAKE_CHECK'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','If the following features are found in coordinate file then make restraints to maintain them:'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'widget','MAKE_PEPTIDE',
        'label','D-peptide ',
        'widget','MAKE_CISPEPTIDE',
        'label','cis-peptide',
        'widget','MAKE_SYMMETRY',
        'label','Links between symmetry related atoms'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','Make links between:'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','Amino acids and DNA/RNA if',
        'widget','MAKE_CONNECTIVITY'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','Sugar-sugar and sugar-peptide if',
        'widget','MAKE_SUGAR'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','All others if',
        'widget','MAKE_LINK'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------
  def drawSetupGeometricRestraintsReview(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['REVIEW']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['REVIEW']])

    self.createLine([
        'widget','SETUP_GEOMETRIC_RESTRAINTS',
        'advice','Setup Geometric Restraints'])

    self.createLine([
        'label','Checking against dictionary:',
        'widget','REVIEW_MAKE_CHECK'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','If the following features are found in coordinate file then make restraints to maintain them:'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'widget','REVIEW_MAKE_PEPTIDE',
        'label','D-peptide ',
        'widget','REVIEW_MAKE_CISPEPTIDE',
        'label','cis-peptide',
        'widget','REVIEW_MAKE_SYMMETRY',
        'label','Links between symmetry related atoms'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','Make links between:'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','Amino acids and DNA/RNA if',
        'widget','REVIEW_MAKE_CONNECTIVITY'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','Sugar-sugar and sugar-peptide if',
        'widget','REVIEW_MAKE_SUGAR'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.createLine([
        'label','All others if',
        'widget','REVIEW_MAKE_LINK'],
        toggle=['SETUP_GEOMETRIC_RESTRAINTS','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------

# def show_NCS_GROUP_LIST(self):
#   par = self.container.controlParameters
#   gui = self.container.guiParameters
#   return gui.SETUP_NCS_RESTRAINTS and not par.IFAUTONCS._value

# def drawSetupNCSRestraints(self):

#   self.openSubFrame(
#       frame=True,
#       toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])
#   self.closeSubFrame()

#   self.openSubFrame(
#       frame=True,
#       toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])

#   self.createLine([
#       'widget','SETUP_NCS_RESTRAINTS',
#       'advice','Setup Non-Crystallographic Symmetry (NCS) Restraints'])

#   self.createLine([
#       'widget','IFAUTONCS',
#       'label','use automatically generated',
#       'widget','AUTONCS_MODE',
#       'label','NCS restraints'],
#       toggle=['SETUP_NCS_RESTRAINTS','open',[True]])

#   self.createLine([
#       'widget','NCS_GROUP_LIST'],
#       toggleFunction=[self.show_NCS_GROUP_LIST,['SETUP_NCS_RESTRAINTS','IFAUTONCS']])

#   self.closeSubFrame()

# ----------------------------------------------------------------------

  def drawSetupNCSRestraints(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])

    self.createLine([
        'widget','SETUP_NCS_RESTRAINTS',
        'advice','Setup Non-Crystallographic Symmetry (NCS) Restraints'])

    self.createLine([
        'widget','IFAUTONCS',
        'label','use automatically generated',
        'widget','AUTONCS_MODE',
        'label','NCS restraints'],
        toggle=['SETUP_NCS_RESTRAINTS','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def show_EXTERNAL_NCYCLES(self):
    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.REFINEMENT_PARAMETERS and par.RUN_COOT_FW._value

  def show_MATRIX_WEIGHT(self):
    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.REFINEMENT_PARAMETERS and not par.AUTO_WEIGHTING._value

  def show_PHASE_BBLUR(self):
    inp = self.container.inputData
    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.REFINEMENT_PARAMETERS and str(par.INPUT_PHASE) in ['PHASE','HL']

# ----------------------------------------------------------------------

  def drawRefinementParametersREST(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['REST']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['REST']])

    self.createLine([
        'widget','REFINEMENT_PARAMETERS',
        'advice','Refinement Parameters'])

    self.createLine([
        'label','Do',
        'widget','NCYCLES',
        'label','cycles of maximum likelihood restrained refinement'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','   in each Refmac run and',
        'widget','EXTERNAL_NCYCLES',
        'label','cycle(s) of Coot:findwaters'],
        toggleFunction=[self.show_EXTERNAL_NCYCLES,['REFINEMENT_PARAMETERS','RUN_COOT_FW']])

    self.createLine([
        'label','Use hydrogen atoms:',
        'widget','MAKE_HYDROGEN',
        'label','and',
        'widget','MAKE_HOUT',
        'label','output to coordinate file'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_RESOLUTION',
        'label','Resolution range from minimum',
        'widget','EXCLUDE_RESOLUTION_MIN',
        'label',' to ',
        'widget','EXCLUDE_RESOLUTION_MAX'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','AUTO_WEIGHTING',
        'label','Use automatic weighting',
        'widget','EXPERIMENTAL_WEIGHTING',
        'label','Use experimental sigmas to weight Xray terms'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Use weighting term',
        'widget','MATRIX_WEIGHT'],
        toggleFunction=[self.show_MATRIX_WEIGHT,['REFINEMENT_PARAMETERS','AUTO_WEIGHTING']])

    self.createLine([
        'widget','IFJELLY',
        'label','use jelly-body refinement with sigma',
        'widget','JELLY_SIGMA'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Refine',
        'widget','B_REFINEMENT_MODE',
        'label','temperature factors'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_FREER',
        'label','Exclude data with freeR label',
        'widget','FREERFLAG',
        'label','with value of',
        'widget','EXCLUDE_FREER_VALUE'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Blurring factors for input phase FOM SC',
        'widget','PHASE_SCBLUR',
        'label','B',
        'widget','PHASE_BBLUR'],
        toggleFunction=[self.show_PHASE_BBLUR,['REFINEMENT_PARAMETERS','INPUT_PHASE']])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def drawRefinementParametersUNRE(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['UNRE']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['UNRE']])

    self.createLine([
        'widget','REFINEMENT_PARAMETERS',
        'advice','Refinement Parameters'])

    self.createLine([
        'label','Do',
        'widget','NCYCLES',
        'label','cycles of maximum likelihood unrestrained refinement'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','   in each Refmac run and',
        'widget','EXTERNAL_NCYCLES',
        'label','cycle(s) of Coot:findwaters'],
        toggleFunction=[self.show_EXTERNAL_NCYCLES,['REFINEMENT_PARAMETERS','RUN_COOT_FW']])

    self.createLine([
        'label','Use hydrogen atoms:',
        'widget','MAKE_HYDROGEN',
        'label','and',
        'widget','MAKE_HOUT',
        'label','output to coordinate file'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_RESOLUTION',
        'label','Resolution range from minimum',
        'widget','EXCLUDE_RESOLUTION_MIN',
        'label',' to ',
        'widget','EXCLUDE_RESOLUTION_MAX'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','AUTO_WEIGHTING',
        'label','Use automatic weighting',
        'widget','EXPERIMENTAL_WEIGHTING',
        'label','Use experimental sigmas to weight Xray terms'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Use weighting term',
        'widget','MATRIX_WEIGHT'],
        toggleFunction=[self.show_MATRIX_WEIGHT,['REFINEMENT_PARAMETERS','AUTO_WEIGHTING']])

    self.createLine([
        'label','Refine',
        'widget','B_REFINEMENT_MODE',
        'label','temperature factors'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_FREER',
        'label','Exclude data with freeR label',
        'widget','FREERFLAG',
        'label','with value of',
        'widget','EXCLUDE_FREER_VALUE'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Blurring factors for input phase FOM SC',
        'widget','PHASE_SCBLUR',
        'label','B',
        'widget','PHASE_BBLUR'],
        toggleFunction=[self.show_PHASE_BBLUR,['REFINEMENT_PARAMETERS','INPUT_PHASE']])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def drawRefinementParametersRIGID(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['RIGID']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['RIGID']])

    self.createLine([
        'widget','REFINEMENT_PARAMETERS',
        'advice','Refinement Parameters'])

    self.createLine([
        'label','Do',
        'widget','RIGID_NCYCLES',
        'label','cycles of maximum likelihood rigid body refinement'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Use hydrogen atoms:',
        'widget','MAKE_HYDROGEN',
        'label','and',
        'widget','MAKE_HOUT',
        'label','output to coordinate file'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_RESOLUTION',
        'label','Resolution range from minimum',
        'widget','EXCLUDE_RESOLUTION_MIN',
        'label',' to ',
        'widget','EXCLUDE_RESOLUTION_MAX'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','AUTO_WEIGHTING',
        'label','Use automatic weighting',
        'widget','EXPERIMENTAL_WEIGHTING',
        'label','Use experimental sigmas to weight Xray terms'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Use weighting term',
        'widget','MATRIX_WEIGHT'],
        toggleFunction=[self.show_MATRIX_WEIGHT,['REFINEMENT_PARAMETERS','AUTO_WEIGHTING']])

    self.createLine([
        'label','Refine overall B-factor'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_FREER',
        'label','Exclude data with freeR label',
        'widget','FREERFLAG',
        'label','with value of',
        'widget','EXCLUDE_FREER_VALUE'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Blurring factors for input phase FOM SC',
        'widget','PHASE_SCBLUR',
        'label','B',
        'widget','PHASE_BBLUR'],
        toggleFunction=[self.show_PHASE_BBLUR,['REFINEMENT_PARAMETERS','INPUT_PHASE']])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def drawRefinementParametersTLS(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['TLS']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['TLS']])

    self.createLine([
        'widget','REFINEMENT_PARAMETERS',
        'advice','Refinement Parameters'])

    self.createLine([
        'label','Do',
        'widget','NCYCLES',
        'label','cycles of maximum likelihood restrained refinement after TLS refinement'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Use hydrogen atoms:',
        'widget','MAKE_HYDROGEN',
        'label','and',
        'widget','MAKE_HOUT',
        'label','output to coordinate file'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_RESOLUTION',
        'label','Resolution range from minimum',
        'widget','EXCLUDE_RESOLUTION_MIN',
        'label',' to ',
        'widget','EXCLUDE_RESOLUTION_MAX'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','AUTO_WEIGHTING',
        'label','Use automatic weighting',
        'widget','EXPERIMENTAL_WEIGHTING',
        'label','Use experimental sigmas to weight Xray terms'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Use weighting term',
        'widget','MATRIX_WEIGHT'],
        toggleFunction=[self.show_MATRIX_WEIGHT,['REFINEMENT_PARAMETERS','AUTO_WEIGHTING']])

    self.createLine([
        'widget','IFJELLY',
        'label','use jelly-body refinement with sigma',
        'widget','JELLY_SIGMA'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Refine',
        'widget','B_REFINEMENT_MODE',
        'label','temperature factors'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'widget','EXCLUDE_FREER',
        'label','Exclude data with freeR label',
        'widget','FREERFLAG',
        'label','with value of',
        'widget','EXCLUDE_FREER_VALUE'],
        toggle=['REFINEMENT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Blurring factors for input phase FOM SC',
        'widget','PHASE_SCBLUR',
        'label','B',
        'widget','PHASE_BBLUR'],
        toggleFunction=[self.show_PHASE_BBLUR,['REFINEMENT_PARAMETERS','INPUT_PHASE']])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def show_PROSMART_STRAND(self):

    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.EXTERNAL_RESTRAINTS and str(par.PROSMART_MODE) == 'SECSTR'

  def drawExternalRestraints(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])

    self.createLine([
        'widget','EXTERNAL_RESTRAINTS',
        'advice','External Restraints'])

    self.createLine([
        'widget','PROSMART_HELIX',
        'label','Apply Prosmart helix restraints, and',
        'widget','PROSMART_STRAND',
        'label','Apply Prosmart strand restraints'],
        toggleFunction=[self.show_PROSMART_STRAND,['EXTERNAL_RESTRAINTS','PROSMART_MODE']])

    self.createLine([
        'widget','IFEXTREST_SCALE',
        'label','Apply external restraints with weight',
        'widget','EXTREST_SCALE',
        'label','and',
        'widget','IFEXTREST_GMWT',
        'label','apply Geman-McClure weight',
        'widget','EXTREST_GMWT'],
        toggle=['EXTERNAL_RESTRAINTS','open',[True]])

    self.createLine([
        'widget','IFEXTREST_DMAX',
        'label','Apply maximum external restraint distance',
        'widget','EXTREST_DMAX',
        'label','and',
        'widget','IFEXTREST_USEMAIN',
        'label',' apply to main chain only'],
        toggle=['EXTERNAL_RESTRAINTS','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------
  def drawIdealisationParameters(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['IDEA']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['IDEA']])

    self.createLine([
        'widget','IDEALISATION_PARAMETERS',
        'advice','Idealisation Parameters'])

    self.createLine([
        'label','Do',
        'widget','NCYCLES',
        'label','cycles of idealisation'],
        toggle=['IDEALISATION_PARAMETERS','open',[True]])

    self.createLine([
        'label','Use hydrogen atoms:',
        'widget','MAKE_HYDROGEN',
        'label','and',
        'widget','MAKE_HOUT',
        'label','output to coordinate file'],
        toggle=['IDEALISATION_PARAMETERS','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------
  def drawCootParameters(self):

    self.openSubFrame(
        frame=True,
        toggle=['RUN_COOT_FW','open',[True]])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['RUN_COOT_FW','open',[True]])

    self.createLine([
        'widget','COOT_PARAMETERS',
        'advice','Coot Parameters'])

    self.createLine([
        'label','Add waters where DELFWT map greater than ',
        'widget','COOT_SIGMA_ADD',
        'label','sigma'],
        toggle=['COOT_PARAMETERS','open',[True]])

    self.createLine([
        'label','Remove waters where FWT map less than ',
        'widget','COOT_SIGMA_REMOVE',
        'label','sigma'],
        toggle=['COOT_PARAMETERS','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------
  def drawTLSParameters(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['TLS']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['TLS']])

    self.createLine([
        'widget','TLS_PARAMETERS',
        'advice','TLS Parameters'])

    self.createLine([
        'label','Number of cycles of TLS refinement',
        'widget','TLS_NCYCLES'],
        toggle=['TLS_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFBFAC_SET',
        'label','Set initial Bfactors to',
        'widget','BFAC_SET',
        'label','(numeric value unimportant)'],
        toggle=['TLS_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFADDU',
        'label','Add TLS contribution to XYZOUT (B factors and ANISOU lines)'],
        toggle=['TLS_PARAMETERS','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def show_MONI_LEVEL(self):

    inp = self.container.inputData
    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.MONITORING_AND_OUTPUT_OPTIONS and not str(par.REFINE_TYPE) == 'RIGID'

  def show_MONI_TORSION(self):

    inp = self.container.inputData
    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.MONITORING_AND_OUTPUT_OPTIONS and not str(par.REFINE_TYPE) == 'UNRE'

  def show_EXTEND_MAP(self):

    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.MONITORING_AND_OUTPUT_OPTIONS and par.IF_MAPOUT._value

# ----------------------------------------------------------------------

  def drawMonitoringAndOutputOptions(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['REVIEW','IDEA']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['REVIEW','IDEA']])

    self.createLine([
        'widget','MONITORING_AND_OUTPUT_OPTIONS',
        'advice','Monitoring and Output Options'])

    self.createLine([
        'widget','IFSHARP',
        'label','Perform map sharpening. Manually set B value',
        'widget','B_SHARP',
        'label','and/or alpha value',
        'widget','AL_SHARP'],
        toggle=['MONITORING_AND_OUTPUT_OPTIONS','open',[True]])

    self.createLine([
        'label','Output',
        'widget','MONI_LEVEL',
        'label','monitoring statistics'],
        toggleFunction=[self.show_MONI_LEVEL,['MONITORING_AND_OUTPUT_OPTIONS','REFINE_TYPE']])

    self.createLine([
        'label','Sigma cutoffs for monitoring levels..'],
        toggleFunction=[self.show_MONI_TORSION,['MONITORING_AND_OUTPUT_OPTIONS','REFINE_TYPE']])

    self.createLine([
        'label','Torsion',
        'label','Distance',
        'label','Angle',
        'label','Planarity',
        'label','vanDerWaals'],
        toggleFunction=[self.show_MONI_TORSION,['MONITORING_AND_OUTPUT_OPTIONS','REFINE_TYPE']])

    self.createLine([
        'widget','MONI_TORSION',
        'widget','MONI_DISTANCE',
        'widget','MONI_ANGLE',
        'widget','MONI_PLANE',
        'widget','MONI_VANDERWAALS'],
        toggleFunction=[self.show_MONI_TORSION,['MONITORING_AND_OUTPUT_OPTIONS','REFINE_TYPE']])

    self.createLine([
        'label','Chiral',
        'label','Bfactor',
        'label','Bsphere',
        'label','Rbond',
        'label','NCSr'],
        toggleFunction=[self.show_MONI_TORSION,['MONITORING_AND_OUTPUT_OPTIONS','REFINE_TYPE']])

    self.createLine([
        'widget','MONI_CHIRAL',
        'widget','MONI_BFACTOR',
        'widget','MONI_BSPHERE',
        'widget','MONI_RBOND',
        'widget','MONI_NCSR'],
        toggleFunction=[self.show_MONI_TORSION,['MONITORING_AND_OUTPUT_OPTIONS','REFINE_TYPE']])

    self.createLine([
        'widget','IF_MAPOUT',
        'label','Generate weighted difference maps files in',
        'widget','MAPOUT_FORMAT',
        'label','format'],
        toggle=['MONITORING_AND_OUTPUT_OPTIONS','open',[True]])

    self.createLine([
        'widget','EXTEND_MAP',
        'label','Extend map to cover molecule with border',
        'widget','MAP_BORDER'],
        toggleFunction=[self.show_EXTEND_MAP,['MONITORING_AND_OUTPUT_OPTIONS','IF_MAPOUT']])

    self.createLine([
        'label','Fwt map',
        'widget','MAPOUT1'],
        toggleFunction=[self.show_EXTEND_MAP,['MONITORING_AND_OUTPUT_OPTIONS','IF_MAPOUT']])

    self.createLine([
        'label','DelFwt map',
        'widget','MAPOUT2'],
        toggleFunction=[self.show_EXTEND_MAP,['MONITORING_AND_OUTPUT_OPTIONS','IF_MAPOUT']])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def show_BULK_SCALING_RESOLUTION_MIN(self):

    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.SCALING and str(par.BULK_SOLVENT_SCALING) == 'BULK'

  def show_SIMPLE_SCALING_RESOLUTION_MIN(self):

    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.SCALING and str(par.BULK_SOLVENT_SCALING) == 'SIMPLE'

  def show_SOLVENT_VDWPROB(self):

    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.SCALING and par.IF_SOLVENT._value

# ----------------------------------------------------------------------

  def drawScaling(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['REVIEW','IDEA']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['REVIEW','IDEA']])

    self.createLine([
        'widget','SCALING',
        'advice','Scaling'])

    self.createLine([
        'label','Use',
        'widget','BULK_SOLVENT_SCALING',
        'label','scaling'],
        toggle=['SCALING','open',[True]])

    self.createLine([
        'label','Bulk solvent scaling for resolution range',
        'widget','BULK_SCALING_RESOLUTION_MIN',
        'label',' to ',
        'widget','BULK_SCALING_RESOLUTION_MAX'],
        toggleFunction=[self.show_BULK_SCALING_RESOLUTION_MIN,['SCALING','BULK_SOLVENT_SCALING']])

    self.createLine([
        'label','Simple solvent scaling for resolution range',
        'widget','SIMPLE_SCALING_RESOLUTION_MIN',
        'label',' to ',
        'widget','SIMPLE_SCALING_RESOLUTION_MAX'],
        toggleFunction=[self.show_SIMPLE_SCALING_RESOLUTION_MIN,['SCALING','BULK_SOLVENT_SCALING']])

    self.createLine([
        'label','Determine scaling using the',
        'widget','SCALING_REF_SET',
        'label','set of reflections   ',
        'widget','SCALING_EXPE_SIGMA',
        'label','Use experimental sigmas'],
        toggle=['SCALING','open',[True]])

    self.createLine([
        'widget','IF_SOLVENT',
        'label','Calculate the contribution from the solvent region'],
        toggle=['SCALING','open',[True]])

    self.createLine([
        'label','For the solvent mask calculation:'],
        toggleFunction=[self.show_SOLVENT_VDWPROB,['SCALING','IF_SOLVENT']])

    self.createLine([
        'label','     Increase VDW radius of non-ion atoms by',
        'widget','SOLVENT_VDWPROB'],
        toggleFunction=[self.show_SOLVENT_VDWPROB,['SCALING','IF_SOLVENT']])

    self.createLine([
        'label','     Increase ionic radius of potential ion atoms by',
        'widget','SOLVENT_IONPROB'],
        toggleFunction=[self.show_SOLVENT_VDWPROB,['SCALING','IF_SOLVENT']])

    self.createLine([
        'label','     Shrink the area of the mask by',
        'widget','SOLVENT_RSHRINK',
        'label','after calculation'],
        toggleFunction=[self.show_SOLVENT_VDWPROB,['SCALING','IF_SOLVENT']])

    self.createLine([
        'widget','SCALING_IF_FIXB',
        'label','For low resolution structures, fix the B values of Babinet bulk solvent to',
        'widget','SCALING_FIXB_BBULK'],
        toggleFunction=[self.show_BULK_SCALING_RESOLUTION_MIN,['SCALING','BULK_SOLVENT_SCALING']])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def drawRigidDomainsDefinition(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['RIGID']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','open',['RIGID']])

    self.createLine([
        'widget','RIGID_DOMAINS_DEFINITION',
        'advice','Rigid Domains Definition'])

    self.createLine([
        'widget','RIGID_GROUP_LIST'],
        toggle=['RIGID_DOMAINS_DEFINITION','open',[True]])

    self.closeSubFrame()

# ----------------------------------------------------------------------

  def show_IFISO(self):

    inp = self.container.inputData
    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.GEOMETRIC_PARAMETERS and str(par.REFINE_TYPE) in ['REST','UNRE'] and str(par.B_REFINEMENT_MODE) in ['ANIS','MIXED']

  def show_BLIM(self):

    inp = self.container.inputData
    par = self.container.controlParameters
    gui = self.container.guiParameters
    return gui.GEOMETRIC_PARAMETERS and not str(par.REFINE_TYPE) == 'IDEA'

# ----------------------------------------------------------------------

  def drawGeometricParameters(self):

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])
    self.closeSubFrame()

    self.openSubFrame(
        frame=True,
        toggle=['REFINE_TYPE','hide',['RIGID','UNRE']])

    self.createLine([
        'widget','GEOMETRIC_PARAMETERS',
        'advice','Geometric parameters'])

    self.createLine([
        'label','Restraint',
        'label','Overall wt',
        'label','Sigmas'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFDIST',
        'label','Distance',
        'widget','WDSKAL'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFANGL',
        'label','Angle',
        'widget','ANGLE_SCALE'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label',' ',
        'label',' ',
        'label','main chain',
        'label','main chain',
        'label','side chain',
        'label','side chain'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label',' ',
        'label',' ',
        'label','bond',
        'label','angle',
        'label','bond',
        'label','angle'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFTMP',
        'label','Bfactor',
        'widget','WBSKAL',
        'widget','SIGB1',
        'widget','SIGB2',
        'widget','SIGB3',
        'widget','SIGB4'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label',' ',
        'label',' ',
        'label','peptide',
        'label','aromatic'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFPLAN',
        'label','Plane',
        'widget','WPSKAL',
        'widget','SIGPP',
        'widget','SIGPA'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label',' ',
        'label',' ',
        'label','chiral'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFCHIR',
        'label','Chiral',
        'widget','WCSKAL',
        'widget','SIGC'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFTORS',
        'label','Torsion',
        'widget','WTSKAL'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label',' ',
        'label',' ',
        'label','tight',
        'label','medium',
        'label','loose'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFNCSR',
        'label','NCS position',
        'widget','WSSKAL',
        'widget','SIGSP1',
        'widget','SIGSP2',
        'widget','SIGSP3'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label','     NCS Bfactor ',
        'label',' ',
        'widget','SIGSB1',
        'widget','SIGSB2',
        'widget','SIGSB3'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label','VDW contacts'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFVAND',
        'label','VDW contacts',
        'widget','WVSKAL'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label','Sigmas for types of non-bonding interaction..'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label','Non-bonding',
        'widget','WAND_SIGMA_VDW',
        'label','H-bonding',
        'widget','WAND_SIGMA_HBOND',
        'label','metal-ion interactions',
        'widget','WAND_SIGMA_METAL',
        'label','1-4 atoms in torsion',
        'widget','WAND_SIGMA_TORS'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label','Set increments for non-bonded interactions..'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'label','1-4 atoms in torsion',
        'widget','WAND_INCR_TORS',
        'label','H-bond pair (not hydrogen atom)',
        'widget','WAND_INCR_ADHB',
        'label','H-bonded pair (one is hydrogen atom)',
        'widget','WAND_INCR_AHHB'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','IFISO',
        'label','Anisotropic refinement',
        'label',' ',
        'label','sphericity',
        'label','bond projection'],
        toggleFunction=[self.show_IFISO,['GEOMETRIC_PARAMETERS','REFINE_TYPE','B_REFINEMENT_MODE']])

    self.createLine([
        'label',' ',
        'label',' ',
        'widget','SPHERICITY',
        'widget','RBOND'],
        toggleFunction=[self.show_IFISO,['GEOMETRIC_PARAMETERS','REFINE_TYPE','B_REFINEMENT_MODE']])

    self.createLine([
        'label','Set limits for B values..'],
        toggle=['GEOMETRIC_PARAMETERS','open',[True]])

    self.createLine([
        'widget','BLIM',
        'label','Limit B value range from',
        'widget','BLIM_MIN',
        'label','to',
        'widget','BLIM_MAX'],
        toggleFunction=[self.show_BLIM,['GEOMETRIC_PARAMETERS','REFINE_TYPE']])

    self.closeSubFrame()

