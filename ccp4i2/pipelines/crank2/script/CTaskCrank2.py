"""
Copyright (C) 2011 University of York, Leiden University
"""

import functools
import os
import traceback

from PySide2 import QtCore

from . import crank2_basepipe
from ....core import CCP4ErrorHandling
from ....core import CCP4Utils
from ....core.CCP4Modules import PROJECTSMANAGER
from ....qtgui import CCP4TaskWidget
from ....qtgui import CCP4Widgets


class CTaskCrank2(CCP4TaskWidget.CTaskWidget):
  TASKNAME = 'crank2'
  TASKVERSION = 0.02
  TASKMODULE='expt_phasing'
  TASKTITLE='Automated structure solution - CRANK2 phasing and building'
  SHORTTASKTITLE='CRANK2'
  DESCRIPTION='CRANK2 experimental phasing pipeline'
  WHATNEXT = ['coot_rebuild','prosmart_refmac',['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml']]
  MGDISPLAYFILES = ['XYZOUT']
  RANK=1

  EXPORTMTZPARAMS = [ ['F_SIGFanom', 'F_SIGFanom2', 'F_SIGFanom3', 'F_SIGFanom4', 'F_SIGFnative'],
                      'FPHOUT_HL', 'FPHOUT_DIFF', 'FPHOUT_2FOFC', 'FPHOUT_DIFFANOM' ]

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)
    # definition of pipelines.
    self.basepipe = crank2_basepipe.crank2_basepipe()
    # these names will appear in the start/end pipeline drop-down menu
    self.step_names = { "substrdet": "Substructure detection", "refatompick": "Substruct. improvement",
                        "phas": "Substructure phasing", "phdmmb": "Den.mod. & poly-Ala tracing",
                        "handdet": "Hand determination", "dmfull": "Density modification", 
                        "building": "Model building", "ref": "Refinement" }

    # a lot of params will need to be added here - params from gui1 can be used as a good starting point!
    self.default_generation_triggers = [ 'SEQIN', 'F_SIGFanom','F_SIGFanom2','F_SIGFanom3', 'F_SIGFanom4','F_SIGFnative', 'NATIVE', 'MAD2',\
                                         'F_SIGFanom_nonmtz','F_SIGFanom2_nonmtz','F_SIGFanom3_nonmtz', 'F_SIGFanom4_nonmtz','F_SIGFnative_nonmtz', \
                                         'USE_COMB', 'SUBSTRDET_PROGRAM', 'SHELXCDE', 'DNAME','DNAME2','DNAME3','DNAME4']
    self.default_generation_string_triggers = [ 'ATOM_TYPE', 'MONOMERS_ASYM', 'WAVELENGTH', 'RESIDUES_MON', \
           'WAVELENGTH2','WAVELENGTH3','WAVELENGTH4','CELL_A','CELL_B','CELL_C','CELL_D','CELL_E','CELL_F','SPACEGROUP' ]
    # a list of parameters that do not obey the crank2 naming style
    # (this list is only used to check for possible typos - not critical, functionality not affected at all)
    self.crank2_naming_exceptions = ("HANDDET_DMCYC","REF_EXCLUDE_FREE","COMB_PHDMMB_EXCLUDE_FREE","MBREF_EXCLUDE_FREE")
    self.save_def_connect, self.save_user_connect = {}, {}
    self.is_def_connect, self.is_user_connect = {}, {}
    self.save_connect = { self.setDefaultParameters: (self.save_def_connect, self.is_def_connect),
                          self.setByUser:            (self.save_user_connect, self.is_user_connect),
                        }
    self.defstr_act = {}
    # saving the input data file path + type to check whether they changed and only generate defaults if they did 
    # to fight with i2 often changing the objects unneceserilly and slowing down the interface responsiveness
    self.prev_data={'F_SIGFanom':(None,None),'F_SIGFanom2':(None,None),'F_SIGFanom3':(None,None),'F_SIGFanom4':(None,None),'F_SIGFnative':(None,None)}
    self.job_started=False


  def getcont(self,strng,cont=None):
    if cont is None:
      cont = self.container
    if self.has_cont_attr(cont.inputData,strng):
      return getattr(cont.inputData,strng)
    elif self.has_cont_attr(cont.controlParameters,strng):
      return getattr(cont.controlParameters,strng)
    else:
      return None

  def has_cont_attr(self,cont,strng):
    try:
      return hasattr(cont,strng)
    except:
      return False

  def ToggleDataInputMtz(self,nonmtz,mad):
    return (mad=='' or self.getcont(mad)) and not self.getcont(nonmtz)

  def ToggleDataInputNonMtz(self,nonmtz,mad):
    return (mad=='' or self.getcont(mad)) and self.getcont(nonmtz)

  def AddInputDataLine(self,fsigf,indent2,nonanom=False):
    si = fsigf[-1]  if fsigf.endswith(('2','3','4'))  else ''
    mad = 'MAD'+si  if si  else ''
    if nonanom:  mad = 'NATIVE'
    if self.TASKNAME!='shelx' or nonanom:
      self.createLine( [ 'label',indent2,'widget',fsigf+'_nonmtz' ], toggleFunction = [functools.partial(self.ToggleDataInputNonMtz,'NON_MTZ',mad), ['NON_MTZ',mad]] )
      self.createLine( [ 'label',indent2,'widget',fsigf ], toggleFunction = [functools.partial(self.ToggleDataInputMtz,'NON_MTZ',mad), ['NON_MTZ',mad]] )
      if not nonanom:
        self.createLine( ['label', indent2+'Anomalous scatter. coef: ','label', "f':",'widget','FPRIME'+si, 'label', 'f":','widget','FDPRIME'+si,'label','wavelength:','widget','WAVELENGTH'+si,'widget','DNAME'+si], toggle = [ mad, 'open', [True]])
    else:
      self.createLine(  [ 'label',indent2,'widget',fsigf+'_nonmtz','widget','DNAME'+si ], toggleFunction = [functools.partial(self.ToggleDataInputNonMtz,'NON_MTZ',mad), ['NON_MTZ',mad]] )
      self.createLine(  [ 'label',indent2,'widget',fsigf,'widget','DNAME'+si ], toggleFunction = [functools.partial(self.ToggleDataInputMtz,'NON_MTZ',mad), ['NON_MTZ',mad]])



  def drawContents(self):
    ctrl, inp = self.container.controlParameters, self.container.inputData
    indent1 = '&nbsp;&nbsp;&nbsp;&nbsp;'
    indent2 = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

    # parameters for which set-by-user is distinguished from crank2 defaults 
    self.set_by_user_params = [ (d[5:],ctrl)  for d in ctrl._dataOrder if d.startswith('USER_') ] + \
                              [ (d[5:],inp)   for d in inp._dataOrder  if d.startswith('USER_') ]

    self.setProgramHelpFile('crank2')

    folder = self.openFolder(folderFunction='inputData',title='Input Data')

    self.createLine(  ['label', 'Start pipeline with', 'widget', 'START_PIPELINE', 'label', 'and end with', 'widget', 'END_PIPELINE'] )
    if self.TASKNAME!='shelx':
      self.createLine(  ['subtitle','Input partial model (MR-SAD, model rebuilding)','If a partial protein model is available from molecular replacement or preliminary model building (for SAD only)','widget','INPUT_PARTIAL'] )
      self.openSubFrame( toggle = [ 'INPUT_PARTIAL', 'open', [True]] )
      self.createLine(  [  'label', indent1,'widget','XYZIN' ] )
      self.createLine(  [  'label', indent1+'Start from anomalous substructure only - remove all non-substructure atoms','widget','PARTIAL_AS_SUBSTR' ] )
      self.closeSubFrame()
      #self.openSubFrame( toggle = [ 'INPUT_PARTIAL', 'open', [False]] )
      #self.closeSubFrame()

    self.createLine(  [ 'subtitle', 'Input protein sequence ', 'widget', 'INPUT_SEQUENCE' ] )
    self.openSubFrame( toggle = [ 'INPUT_SEQUENCE', 'open', [True]] )
    self.createLine(  [ 'label',indent1,'widget','SEQIN' ] )
    infoline=self.createLine(  [ 'advice', indent1+indent1+'<small>In total</small>', ], toggleFunction = [self.ToggleSeq, ['SEQIN','RESIDUES_MON']] )
    infoline.addWidget( CCP4Widgets.CStringView(self, model=inp.RESIDUES_MON_INFO, qualifiers={'editable':False}) )
    self.createLine(  [ 'advice', '<small>residues found in the sequence.</small>' ], appendLine=infoline )
    self.closeSubFrame()
    self.openSubFrame( toggle = [ 'INPUT_SEQUENCE', 'open', [False]] )
    self.createLine(  ['label', indent1+'Number of residues per monomer', 'widget', 'RESIDUES_MON_COPY'])
    self.closeSubFrame()

    #self.openSubFrame( toggle = [ 'DNA', 'open', [True]] )
    #self.createLine(  ['label', 'Number of nucleotides per monomer', 'widget', 'NUCLEOTIDES_MON'])
    #self.closeSubFrame()

    self.createLine( ['subtitle','Input starting phases','widget','INPUT_PHASES'], toggleFunction = [self.TogglePhases, ['START_PIPELINE','INPUT_PARTIAL']] )
    self.openSubFrame( toggle=['INPUT_PHASES','open',[True]] )
    self.createLine(  ['widget','FPHIN_HL']  )
    self.closeSubFrame()

    self.createLine(['subtitle', 'Crystal #1 composition and collected anomalous dataset(s)', 'Specify the anomalously scattering atoms and dataset(s) for crystal used for phasing.'])
    self.openSubFrame(frame=False)
    
    #at_line=self.createLine(['label', indent1+'Substructure atom:','widget','ATOM_TYPE'])
    #self.createLine( ['label', indent1+'Number of substructure atoms in asym. unit:','widget','NUMBER_SUBSTRUCTURE'], \
    #                 appendLine=at_line, toggle = [ 'ATOM_TYPE', 'close', ['S']] )
    #self.createLine( ['label', ' Sulphurs in asym. unit:','widget','NUMBER_SUBSTRUCTURE',\
    #                  'label',', forming ','widget','SUBSTRDET_NUM_DSUL','label','disulphides'], \
    #                 appendLine=at_line, toggle = [ 'ATOM_TYPE', 'open', ['S']] )
    at_line=self.createLine(['label', indent1+'Substructure atom:','widget','ATOM_TYPE','label','Number of substr. atoms in asymmetric unit:','widget','NUMBER_SUBSTRUCTURE'])
    self.createLine( ['label',indent1+'Number of S-S pairs searched for as 1 supersulfur: ','widget','SUBSTRDET_NUM_DSUL'], toggleFunction = [self.ToggleDSUL, ['ATOM_TYPE','START_PIPELINE']] )
    #at_line.addWidget( CCP4Widgets.CStringView(self, model=inp.ATOM_TYPE, qualifiers={'editable':False}) )
    #self.createLine( ['label','in AU:','widget','NUMBER_SUBSTRUCTURE'], appendLine=at_line )
    #self.createLine( ['label',', forming ','widget','NUMBER_SUBSTRUCTURE_DISULPHIDES','label','disulphides'], appendLine=at_line, toggle = [ 'ATOM_TYPE', 'open', ['S']] )
    self.openSubFrame(toggleFunction = [self.ToggleSubstrModel, ['START_PIPELINE',]])
    self.createLine( ['label',indent1+'Substructure','widget','XYZIN_SUB' ] )
    self.closeSubFrame()
    self.openSubFrame( toggle = [ 'NON_MTZ', 'open', [True]] )
    self.createLine( ['label',indent1+'Cell:','widget','CELL_A','widget','CELL_B','widget','CELL_C',\
                       'widget','CELL_D','widget','CELL_E','widget','CELL_F','label',' Spacegroup:','widget','SPACEGROUP'] )
    self.closeSubFrame()
    self.openSubFrame(frame=True)
    self.createLine( [ 'advice',indent1+'Anomalous data (Friedel pairs)','label',indent1+'Input unmerged/merged SCA/XDS/SHELX format','widget','NON_MTZ' ] )
    self.AddInputDataLine('F_SIGFanom',indent2)
    self.closeSubFrame()
    self.openSubFrame(frame=True, toggle=['INPUT_PARTIAL','close',[True]])
    self.createLine( [ 'label', indent1+'Input anomalous data #2 (MAD)', 'widget', 'MAD2'] )
    self.AddInputDataLine('F_SIGFanom2',indent2)
    self.closeSubFrame()
    self.openSubFrame(frame=True, toggleFunction = [self.ToggleMAD3input, ['INPUT_PARTIAL','MAD2','MAD3']])
    self.createLine( [ 'label', indent1+'Input anomalous data #3 (MAD)', 'widget', 'MAD3'] )
    self.AddInputDataLine('F_SIGFanom3',indent2)
    self.closeSubFrame()
    self.openSubFrame(frame=True, toggleFunction = [self.ToggleMAD4input, ['INPUT_PARTIAL','MAD3','MAD4']])
    self.createLine( [ 'label', indent1+'Input anomalous data #4 (MAD)', 'widget', 'MAD4'] )
    self.AddInputDataLine('F_SIGFanom4',indent2)
    self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Input native observations (Crystal #2)', 'Specify "native" data from a crystal to be modelled, if available', 'widget', 'NATIVE'] )
    self.openSubFrame(frame=True,toggle = ['NATIVE', 'open', [True] ])
    self.AddInputDataLine('F_SIGFnative',indent2='',nonanom=True)
    #self.openSubFrame(toggle = ['NATIVE', 'open', [True] ] )
    self.createLine( [ 'label', indent1+'Substructure atoms present in the native crystal', 'widget', 'SUBSTR_ATOMS_NATIVE' ] )
    self.closeSubFrame()
    freeline = self.createLine( ['subtitle', 'Exclude free ','widget','FREE'], toggleFunction = [self.ToggleFree,['END_PIPELINE']] )
    self.createLine( ['label',' :','widget','FREERFLAG'], appendLine=freeline, toggle = ['FREE', 'open', ['existing'] ] )
    self.createLine( ['label',' constisting of ','widget','FREE_RATIO','label', '% of reflections'], appendLine=freeline, toggle = ['FREE', 'open', ['new'] ] )

    self.openFolder(title='Important Options', drawFolder=self.drawImportant)
    self.openFolder(title='Advanced Options', drawFolder=self.drawAdvanced)

    if self.TASKNAME=='shelx' and not inp.ATOM_TYPE.isSet():
      inp.ATOM_TYPE.set('Se')
    # sending the signal to mark set-by-user
    # qt apparently sends this also when changed manually by defaults and blocking signals did not work,
    # so we reset afterwards in that case at the end of SetI2DefPar
    for p,cont in self.set_by_user_params:
      #getattr(cont,p).dataChanged'.connect(self.setByUser(self,p,cont))
      #getattr(cont,p).dataChanged.connect(functools.partial(self.setByUser,p))
      self.ConnectOneTrig(p,getattr(cont,p),partial=True,funct=self.setByUser)

    # use FREE_EXISTING by default if supplied by i2 (is this safe?) - this certainly cannot be done here as cloned job could start with different assignments!
    #if self.isEditable() and bool(str(inp.FREERFLAG).strip()):
    #  inp.FREE_EXISTING.set(True),  inp.FREE_NEW.set(False)
    # SAD/MAD
    inp.EXPTYPE.dataChanged.connect(self.SetSAD )
    inp.MAD2.dataChanged.connect(self.SADorMAD )
    inp.MAD3.dataChanged.connect(self.SADorMAD )
    inp.MAD4.dataChanged.connect(self.SADorMAD )
    inp.NON_MTZ.dataChanged.connect(self.SADorMAD )
    inp.NON_MTZ.dataChanged.connect(self.CellSpacegroup )
    inp.NATIVE.dataChanged.connect(self.InputNative )
    inp.ATOM_TYPE.dataChanged.connect(self.InputNative )
    inp.NON_MTZ.dataChanged.connect(self.InputNative )
    inp.ATOM_TYPE.dataChanged.connect(self.SSAD )
    inp.START_PIPELINE.dataChanged.connect(self.InputModelPhases )
    inp.START_PIPELINE.dataChanged.connect(self.AdjustEndPipeline )
    inp.END_PIPELINE.dataChanged.connect(self.SetFree )
    inp.SHELXCDE.dataChanged.connect(self.SetPipeline )
    inp.INPUT_PARTIAL.dataChanged.connect(self.SetPipeline )
    inp.INPUT_PARTIAL.dataChanged.connect(self.InputModelPhases )
    inp.EXPTYPE.dataChanged.connect(self.SetPipeline )
    # sequence inputted or not
    inp.INPUT_SEQUENCE.dataChanged.connect(self.InputSequence )
    # the residues_mon from basic options is equal to residues_mon from input data
    # (the reason why they are separated is the fact that i2 would draw the red boxes in basic options but we need them in input data)
    inp.RESIDUES_MON_COPY.dataChanged.connect(self.InputResidues )
    inp.RESIDUES_MON.dataChanged.connect(self.InputResidues2 )

    # connect the defaults generation triggers
    self.ConnectDefaultGenTrig()
    self.ConnectDefaultGenStringTrig()

    # to make sure a new job is properly initialized incl. the red boxes
    if self.isEditable():
      self.AdjustEndPipeline()
      self.CellSpacegroup()
      self.InputSequence(init=True)
      self.SADorMAD(init=True)
      self.InputNative(init=True)
      self.SSAD()
      self.InputModelPhases()

    #folder = self.openFolder(folderFunction='controlParameters',title='Important Options')


  def drawImportant(self):

    ctrl, inp = self.container.controlParameters, self.container.inputData
    indent1 = '&nbsp;&nbsp;&nbsp;&nbsp;'
    indent2 = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

    self.createLine(['label', 'Residues/monomer:','widget', 'RESIDUES_MON', \
                     'label', ' NCS copies:', 'widget', 'MONOMERS_ASYM', \
                     'label', ' Solvent Content: ', 'widget', 'SOLVENT_CONTENT' ])
    self.openSubFrame(toggle = ['DNA', 'open', [True] ] )
    self.createLine(['label', 'Nucleotides per monomer:','widget', 'NUCLEOTIDES_MON'])
    self.closeSubFrame()
    #self.createLine( ['widget','REPLACE_MET_MSE','label','Replace methionones by selenomethionine in the input partial model'])
    self.openSubFrame(toggleFunction = [self.ToggleSADSIRASOption,['EXPTYPE','NATIVE']])
    self.setMenuText('EXPTYPE',{'SAD':'SAD','SIRAS':'SIRAS','MAD':''})
    self.createLine(  [ 'label', 'Use ', 'widget', 'EXPTYPE', 'label', ' phasing and phase improvement'] )
    self.closeSubFrame()
    #self.createLine(['label',''])
    #self.createLine(['label','Job output presentation style ','widget','PRESENT_STYLE'])

    #self.createLine(['label', 'exclude free reflection as defined in UNSSIGNED'])
    #self.createLine(['label', 'Create new free set to exclude ', 'label', ' of reflections'])

    # connect the defaults generation triggers - needs to be done again due to separate tabs
    self.ConnectDefaultGenTrig()
    self.ConnectDefaultGenStringTrig()


  def drawAdvanced(self):

    #folder = self.openFolder(folderFunction='controlParameters',title='Advanced Options')
    ctrl, inp = self.container.controlParameters, self.container.inputData
    indent1 = '&nbsp;&nbsp;&nbsp;&nbsp;'
    indent2 = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

    self.createLine(['advice', indent2+'  <i>Note: empty input fields mean that internal program defaults will be used.</i>'])
    if self.TASKNAME!='shelx':
      # probably input partial not true needs to be added here!
      self.createLine(  ['label', 'Use SHELXC/D/E ','widget','SHELXCDE'], toggleFunction = [self.ToggleSHELXCDE,['START_PIPELINE','INPUT_PARTIAL']] )
    substrdet_toggle,substrdet_toggle_not = [self.basepipe.ToggleDetection,['START_PIPELINE','INPUT_PARTIAL']], [self.basepipe.ToggleNotDetection,['START_PIPELINE','INPUT_PARTIAL']]
    self.openSubFrame(frame=False, toggleFunction = substrdet_toggle)
    self.createLine(['subtitle', 'Substructure detection'])
    if self.TASKNAME!='shelx':
      self.createLine(['label', indent2+'FA Estimation program: ', 'widget', 'FAEST_PROGRAM', 'label', indent1+'Detection program: ', 'widget', 'SUBSTRDET_PROGRAM'])
    self.createLine(['label', indent2+'Initial high resol. cutoff:', 'widget', 'SUBSTRDET_HIGH_RES_CUTOFF', 'label', indent1+'Cutoff radius:','widget','SUBSTRDET_HIGH_RES_CUTOFF_RADIUS', 'label', indent1+'Cutoff step:', 'widget', 'SUBSTRDET_HIGH_RES_CUTOFF_STEP'], \
                    toggle = ['SUBSTRDET_PROGRAM', 'open', ['prasa']] )
    self.createLine(['label', indent2+'High resolution cutoff:', 'widget', 'SUBSTRDET_HIGH_RES_CUTOFF'], \
                    toggle = ['SUBSTRDET_PROGRAM', 'open', ['shelxd','crunch2']] )
    self.createLine(['label', indent2+'Use CCanom1/2 based cutoff (if available):', 'widget', 'SUBSTRDET_HIGH_RES_CUTOFF_CCHALF'], \
                    toggleFunction = [self.ToggleCChalf_rescut,['FAEST_PROGRAM','NON_MTZ']] )
    plin=self.createLine(['label', indent2+'Num. trials:','widget','SUBSTRDET_NUM_TRIALS', 'label', indent1+'CFOM threshold:', 'widget', 'SUBSTRDET_THRESHOLD_STOP'] )
    self.createLine( ['label', indent1+'Min.atoms','widget','PRASA_MINPEAKS', 'label', indent1+'Max.atoms','widget','SUBSTRDET_NUM_ATOMS'], \
                     appendLine=plin, toggle = [ 'SUBSTRDET_PROGRAM', 'open', ['prasa']] )
    #self.createLine( ['label', indent1+'','widget','SUBSTRDET_NUM_ATOMS'], \
    #                 appendLine=plin2, toggle = [ 'PRASA_NUM_ATOMS_RESTR', 'open', [True]] )
    opt1_line = self.createLine(['label', indent2+'Minimum distance between atoms: ', 'widget', 'SUBSTRDET_MIN_DIST_ATOMS', 'label', '&#8491;'+indent1+'Atoms in special positions allowed','widget','SUBSTRDET_MIN_DIST_SYMM_ATOMS'], toggle = ['SUBSTRDET_PROGRAM', 'close', ['crunch2'] ])
    self.createLine(['label', indent1+'Optimize solutions','widget','SUBSTRDET_OPTIMIZE_SOL'], toggle = ['SUBSTRDET_PROGRAM', 'open', ['prasa'] ], appendLine=opt1_line)
    #self.createLine(['label', indent2+'CFOM threshold: ', 'widget', 'SUBSTRDET_THRESHOLD_STOP', 'label', ' Number of trials:','widget','SUBSTRDET_NUM_TRIALS'])
    self.createLine( ['label', indent2+'Number of CPU threads','widget','SUBSTRDET_NUM_THREADS'], \
                     toggle = [ 'SUBSTRDET_PROGRAM', 'open', ['shelxd','prasa']] )
    self.createLine(['label', indent2+'Custom program keywords (comma separated)', 'widget', 'KEYWORDS_SUBSTRDET'])
    self.closeSubFrame()

    self.openSubFrame(frame=False, toggleFunction = [self.basepipe.TogglePeakSearch,['START_PIPELINE','INPUT_PARTIAL','EXPTYPE','SHELXCDE']])
    self.createLine(['subtitle', 'Substructure improvement'])
    #self.createLine(['label', indent2+'Refinement program: ', 'widget', 'REFATOMPICK_REF_PROGRAM'])
    #self.createLine(['label', 'Keywords ', 'widget', 'KEYWORDS_PEAK'])
    self.createLine(['label', indent2+'Max. num. of iterations: ', 'widget', 'REFATOMPICK_NUM_ITER', 'label', indent1+'Number of ref. cycles per iteration ','widget','REFATOMPICK_REFCYC'])
    self.createLine(['label', indent2+'Pick new atoms from anom. maps at peaks above RMS: ', 'widget', 'REFATOMPICK_RMS_THRESHOLD', 'label', indent1+'Remove atoms with occupancy below: ', 'widget', 'REFATOMPICK_OCC_CUT'])
    self.closeSubFrame()

    self.openSubFrame(frame=False, toggleFunction=[self.basepipe.TogglePhasing,['END_PIPELINE','SHELXCDE','INPUT_PARTIAL','PARTIAL_AS_SUBSTR','EXPTYPE']])
    self.createLine(['subtitle', 'Substructure phasing'])
    self.createLine(['label', indent2+'Phasing program: ', 'widget', 'PHAS_PROGRAM'])
    #self.createLine(['label', 'Keywords', 'widget', 'KEYWORDS_PHAS'])
    self.closeSubFrame()

    self.openSubFrame(frame=False, toggleFunction=[self.basepipe.ToggleHandDetermination,['END_PIPELINE','SHELXCDE','INPUT_PARTIAL','PARTIAL_AS_SUBSTR','DO_HANDDET']])
    self.createLine(['subtitle', 'Hand determination', 'widget','DO_HANDDET'])
    #self.createLine(['label', indent2+'Number of iterations: ', 'widget', 'HANDDET_DMCYC'])
    self.createLine(['label', indent2+'Requested hand determination discrimination: ', 'widget', 'HANDDET_THRESHOLD_DISCRIM'], \
           toggleFunction=[self.ToggleHandDeterminationDoHand,['END_PIPELINE','SHELXCDE','INPUT_PARTIAL','PARTIAL_AS_SUBSTR','DO_HANDDET']])
    self.createLine(['label', indent2+'Density modif. program:', 'widget', 'HANDDET_DMFULL_DM_PROGRAM','label', indent1+'Phase combination program: ', 'widget', 'HANDDET_DMFULL_PHCOMB_PROGRAM'], \
           toggleFunction=[self.ToggleHandDeterminationDoHand,['END_PIPELINE','SHELXCDE','INPUT_PARTIAL','PARTIAL_AS_SUBSTR','DO_HANDDET']])
    #self.createLine(['label', 'Keywords for dens.mod.program', 'widget', 'KEYWORDS_HANDDET_DM'])
    #self.createLine(['label', indent2+'Phase combination program: ', 'widget', 'HANDDET_DMFULL_PHCOMB_PROGRAM'])
    #   self.createLine(['label', 'Keywords for Phase combination program', 'widget', 'KEYWORDS_HANDDET_PHCOMB'])
    self.closeSubFrame()

    self.openSubFrame(frame=False, toggleFunction=[self.basepipe.ToggleDensityModification,['END_PIPELINE','SHELXCDE','INPUT_PARTIAL','PARTIAL_AS_SUBSTR']])
    self.createLine(['subtitle', 'Density modification'])
    self.createLine(['label', indent2+'Density modif. program:', 'widget', 'DMFULL_DM_PROGRAM','label', indent1+'Phase combination program: ', 'widget', 'DMFULL_PHCOMB_PROGRAM'])
    #self.createLine(['label', indent2+'Phase combination program: ', 'widget', 'DMFULL_PHCOMB_PROGRAM'])
    self.createLine(['label', indent2+'Number of iterations: ', 'widget', 'DMFULL_DMCYC', 'label', indent1+'FOM threshold: ', 'widget', 'DMFULL_THRESHOLD_STOP'])
    self.createLine(['label', indent2+'Custom options for dens.mod. program:', 'widget', 'KEYWORDS_DMFULL_DM'])
    #    self.createLine(['label', 'Keywords for Phase combination program', 'widget', 'KEYWORDS_DMFULL_PHCOMB'])
    self.closeSubFrame()

    self.openSubFrame(frame=False, toggleFunction=[self.basepipe.ToggleShelxCDE,['END_PIPELINE','SHELXCDE','INPUT_PARTIAL']] )
    self.createLine(['subtitle', 'Density modification and poly-Ala tracing with SHELXE'])
    self.createLine(['label', indent2+'Number of density modif. cycles: ', 'widget', 'PHDMMB_DMCYC', \
                     'label', indent1+'Number of model building cycles: ', 'widget', 'PHDMMB_BIGCYC'])
    self.createLine(['label', indent2+'CC threshold: ', 'widget', 'PHDMMB_THRESHOLD_STOP', 'label', '  Other hand CC threshold: ', 'widget', 'PHDMMB_THRESHOLD_HAND_STOP'])
    self.createLine(['label', indent2+'Use thorough building if CFOM from substr. detection is smaller than', 'widget', 'SUBSTRDET_THRESHOLD_WEAK'], toggleFunction=substrdet_toggle)
    self.createLine(['label', indent2+'Use thorough building', 'widget', 'PHDMMB_THOROUGH_BUILD'], toggleFunction=substrdet_toggle_not)
    self.createLine(['label', indent2+'Custom program arguments (comma separated)', 'widget', 'ARGUMENTS_SHELXE'])
    self.closeSubFrame()

    self.openSubFrame(frame=False,toggleFunction=[self.basepipe.ToggleModelBuilding,['END_PIPELINE','INPUT_PARTIAL']])
    self.createLine(['subtitle', 'Model building'])
    self.createLine(['label', indent2+'Combine phase, model and density modif. information ', 'widget', 'USE_COMB'], toggle = [ 'MB_PROGRAM', 'close', ['arpwarp']] )
    mb_prog_line=self.createLine(['label', indent2+'Model building program: ', 'widget', 'MB_PROGRAM'])
    mb_toggle=[self.basepipe.ToggleUseComb,['USE_COMB','END_PIPELINE','INPUT_PARTIAL','MB_PROGRAM']]
    self.createLine(['label', indent2+'Density modif. program: ', 'widget', 'COMB_PHDMMB_DMFULL_DM_PROGRAM'], toggleFunction=mb_toggle, appendLine=mb_prog_line)
    self.createLine(['label', indent2+'Custom keywords for building program:', 'widget', 'KEYWORDS_MB'])
    self.createLine(['label', indent2+'Custom keywords for dens.mod. program:', 'widget', 'KEYWORDS_COMB_DM'], toggleFunction=mb_toggle)
    self.createLine(['label', indent2+'Start with a few SHELXE tracing cycles', 'widget', 'COMB_PHDMMB_START_SHELXE', 'label', 'with custom keywords', 'widget', 'KEYWORDS_COMB_SHELXE'], toggleFunction=mb_toggle)
    #self.createLine(['label', 'Keywords for building program:', 'widget', 'KEYWORDS_MB'])
    #self.createLine(['label', 'Keywords for refinement program:', 'widget', 'KEYWORDS_MB_REF'])
    self.openSubFrame(frame=False, toggleFunction=mb_toggle)
    #self.createLine(['label', indent2+'Refinement program: ', 'widget', 'COMB_PHDMMB_DMFULL_REF_PROGRAM'])
    #self.createLine(['label', indent2+'DM program: ', 'widget', 'COMB_PHDMMB_DMFULL_DM_PROGRAM'])
    self.createLine(['label', indent2+'Minimum number of cycles: ', 'widget', 'COMB_PHDMMB_MINBIGCYC', 'label', indent1+'Maximum number of cycles: ', 'widget', 'COMB_PHDMMB_MAXBIGCYC'])
    ncs_line=self.createLine(['label', indent2+'Try to determine NCS',  'widget', 'COMB_PHDMMB_NCS_DET'], toggleFunction=[self.basepipe.ToggleNCS,['COMB_PHDMMB_DMFULL_DM_PROGRAM','MONOMERS_ASYM'] ] )
    self.createLine(['label', indent1+'from partial model (rather than heavy atoms)',  'widget', 'COMB_PHDMMB_NCS_DET_MR'], \
             appendLine=ncs_line, toggle = [ 'COMB_PHDMMB_NCS_DET', 'open', [True]] )
    self.createLine(['label', indent2+'Parallel model building:  Use', 'widget', 'COMB_PHDMMB_NUM_PARALLEL', 'label','simultanous building and refinement processes'])
    self.createLine(['label', indent2+'Skip the first model building cycle', 'widget', 'COMB_PHDMMB_SKIP_INITIAL_BUILD', 'label',indent1+'Soft rebuilding','widget','COMB_PHDMMB_REBUILD_ONLY'])
    self.createLine(['label', indent2+'Exclude the free reflections in model building ', 'widget', 'COMB_PHDMMB_EXCLUDE_FREE'], toggle = ['FREE', 'close', ['no']] )
    self.closeSubFrame()
    self.openSubFrame(frame=False, toggleFunction=[self.ToggleNoUseComb,['USE_COMB','END_PIPELINE','INPUT_PARTIAL','MB_PROGRAM']])
    #self.createLine(['label', indent2+'Refinement program: ', 'widget', 'MBREF_REF_PROGRAM'])
    self.createLine(['label', indent2+'Number of building cycles: ', 'widget', 'MBREF_BIGCYC'])
    self.createLine(['label', indent2+'Exclude the free reflections in model building ', 'widget', 'MBREF_EXCLUDE_FREE'], toggle = ['FREE', 'close', ['no']] )
    self.closeSubFrame()
    self.closeSubFrame()

    self.openSubFrame(frame=False,toggleFunction=[self.basepipe.ToggleRefine,['END_PIPELINE','INPUT_PARTIAL']])
    self.createLine(['subtitle', 'Final model refinement with REFMAC'])
    #self.createLine(['label', indent2+'Refinement program: ', 'widget', 'REF_PROGRAM'])
    self.createLine(['label', indent2+'Number of refinement cycles: ', 'widget', 'REF_CYCLES'])
    self.createLine(['label', indent2+'Exclude the free reflections in refinement ', 'widget', 'REF_EXCLUDE_FREE'], toggle = ['FREE', 'close', ['no']] )
    #self.createLine(['label', 'Keywords for refinement program:', 'widget', 'KEYWORDS_REF'])
    #self.createLine(['advice', 'Separate any keywords with commas'])
    self.closeSubFrame()
    self.createLine(['label', 'Remove all intermediate mtz files at the end?', 'widget', 'CLEANUP'])

    # connect the defaults generation triggers - needs to be done again due to separate tabs
    self.ConnectDefaultGenTrig()
    self.ConnectDefaultGenStringTrig()



  @QtCore.Slot()
  def InputResidues(self):
    self.container.inputData.RESIDUES_MON = self.container.inputData.RESIDUES_MON_COPY

  @QtCore.Slot()
  def InputResidues2(self):
    self.container.inputData.RESIDUES_MON_INFO.set( '<small>'+str(self.container.inputData.RESIDUES_MON)+'</small>' )

  @QtCore.Slot(bool)
  def SetPipeline(self,base_steps_only=False):
    self.basepipe.SetBaseSteps(self.container)
    if not base_steps_only:
      inp,ctrl=self.container.inputData,self.container.controlParameters
      self.SetI2DefPar(ctrl, 'DO_HANDDET', False) if inp.INPUT_PARTIAL else self.SetI2DefPar(ctrl, 'DO_HANDDET', True)
      prev_start = str(inp.START_PIPELINE)
      inp.START_PIPELINE.setQualifier('enumerators', list( self.basepipe.base_steps ))
      inp.START_PIPELINE.setQualifier('menuText', list( [self.step_names[s] for s in self.basepipe.base_steps] ))
      if prev_start in inp.START_PIPELINE.qualifiers('enumerators'):
        inp.START_PIPELINE.set(prev_start)
      else:
        inp.START_PIPELINE.set(self.basepipe.base_steps[0])
      if self.getWidget('START_PIPELINE'):
        self.getWidget('START_PIPELINE').validate()
        self.getWidget('START_PIPELINE').populateComboBox(model = inp.START_PIPELINE)
        self.getWidget('START_PIPELINE').updateViewFromModel()
      self.AdjustEndPipeline()

  @QtCore.Slot(str)
  def AdjustEndPipeline(self,step=None):
    # this will need to be changed for MR-SAD and SHELX...
    inp=self.container.inputData
    if step is None:
      step = str(inp.START_PIPELINE)
    prev_end = str(inp.END_PIPELINE)
    end_steps = self.basepipe.base_steps[self.basepipe.base_steps.index(step):]
    inp.END_PIPELINE.setQualifier('enumerators', list( end_steps ))
    inp.END_PIPELINE.setQualifier('menuText', list( [self.step_names[s] for s in end_steps] ))
    if prev_end in inp.END_PIPELINE.qualifiers('enumerators') and prev_end!='ref':
      inp.END_PIPELINE.set(prev_end)
    else:
      #inp.END_PIPELINE.set("ref")
      inp.END_PIPELINE.set("building")
    if self.getWidget('END_PIPELINE'):
      self.getWidget('END_PIPELINE').validate()
      self.getWidget('END_PIPELINE').populateComboBox(model = inp.END_PIPELINE)
      self.getWidget('END_PIPELINE').updateViewFromModel()


  @QtCore.Slot()
  def SetFree(self):
    if not self.ToggleFree():
      self.container.inputData.FREE.set('no')

  @QtCore.Slot()
  def InputModelPhases(self):
    if self.TASKNAME!='shelx':
      self.container.inputData.XYZIN.setQualifier( 'allowUndefined', not bool(self.container.inputData.INPUT_PARTIAL) )
      self.getWidget('XYZIN').validate()
      self.container.inputData.ATOM_TYPE.setQualifier( 'allowUndefined', False )
      self.getWidget('ATOM_TYPE').validate()
    self.container.inputData.XYZIN_SUB.setQualifier( 'allowUndefined', not self.ToggleSubstrModel() or bool(self.container.inputData.INPUT_PARTIAL) )
    self.getWidget('XYZIN_SUB').validate()
    if bool(self.container.inputData.INPUT_PHASES) and not self.TogglePhases():
      self.container.inputData.INPUT_PHASES.set(False)
    self.container.inputData.FPHIN_HL.setQualifier( 'allowUndefined', not self.container.inputData.INPUT_PHASES )
    self.getWidget('FPHIN_HL').validate()

  def ToggleCChalf_rescut(self):
    return str(self.container.controlParameters.FAEST_PROGRAM)=='shelxc' and bool(self.container.inputData.NON_MTZ)

  def ToggleSubstrModel(self):
    return str(self.container.inputData.START_PIPELINE) in ['phdmmb','phas','handdet','dmfull','building','ref','refatompick' ]

  def ToggleDSUL(self):
    return self.basepipe.CheckStartEnd('substrdet') and str(self.container.inputData.ATOM_TYPE)=='S'

  def TogglePhases(self):
    # input of phases to shelxe not supported as of now although technically possible
    return str(self.container.inputData.START_PIPELINE) in ['handdet','dmfull','building','ref' ] or bool(self.container.inputData.INPUT_PARTIAL)

  def ToggleSHELXCDE(self):
    return str(self.container.inputData.START_PIPELINE) in ['substrdet','phas','phdmmb','refatompick'] and not bool(self.container.inputData.INPUT_PARTIAL)

  def ToggleSeq(self):
    return bool(str(self.container.inputData.SEQIN).strip())

  def ToggleMAD3input(self):
    inp = self.container.inputData
    return not inp.INPUT_PARTIAL and (inp.MAD2 or inp.MAD3)

  def ToggleMAD4input(self):
    inp = self.container.inputData
    return not inp.INPUT_PARTIAL and (inp.MAD3 or inp.MAD4)

  @QtCore.Slot(bool)
  def InputSequence(self,init=False):
    inp=self.container.inputData
    inp.SEQIN.setQualifier( 'allowUndefined', not bool(inp.INPUT_SEQUENCE) )
    self.getWidget('SEQIN').validate()
    inp.RESIDUES_MON_COPY.setQualifier( 'mustExist', not bool(inp.INPUT_SEQUENCE) )
    inp.RESIDUES_MON_COPY.setQualifier( 'allowUndefined', bool(inp.INPUT_SEQUENCE) )
    self.getWidget('RESIDUES_MON_COPY').validate()
    if not init and inp.INPUT_SEQUENCE:
      inp.USER_RESIDUES_MON.set(False)

  def ConnectDefaultGenStringTrig(self):
    for trig in self.default_generation_string_triggers:
      #self.getWidget(trig.objectName()).widget.editingFinished.connect(self.setDefaultParameters )
      if self.getWidget(trig) and hasattr(self.getWidget(trig),'widget'):
        self.defstr_act[trig] = [False, False]
        self.getWidget(trig).widget.editingFinished.connect(functools.partial(self.DefStrTrig0,trig) )
        self.getcont(trig).dataChanged.connect(functools.partial(self.DefStrTrig1,trig) )

  # making sure that both editingFinished and dataChanged is connected to prevent useless def. generation
  @QtCore.Slot(str)
  def DefStrTrig0(self,trig):
    self.defstr_act[trig][0] = True
    if self.defstr_act[trig][1]:
      self.defstr_act[trig][0], self.defstr_act[trig][1] = False, False
      self.setDefaultParameters(trig)

  @QtCore.Slot(str)
  def DefStrTrig1(self,trig):
    self.defstr_act[trig][1] = True
    if self.defstr_act[trig][0]:
      self.defstr_act[trig][0], self.defstr_act[trig][1] = False, False
      self.setDefaultParameters(trig)

  def ConnectOneTrig(self,trig,obj,disconnect=False,partial=True,funct=None):
    # assuming default generation triggers by default
    if funct is None:
      funct = self.setDefaultParameters
    saved_conn = self.save_connect[funct][0]
    is_conn = self.save_connect[funct][1]
    if partial and trig not in saved_conn:
      saved_conn[trig] = functools.partial(funct,trig)
    defpar = saved_conn[trig]  if partial and trig in saved_conn  else funct
    if disconnect or trig not in is_conn or not is_conn[trig]:
      if disconnect:
        if trig in is_conn:
          obj.dataChanged.disconnect(defpar)
        is_conn[trig] = False
      else:
        is_conn[trig] = True
        obj.dataChanged.connect(defpar)

  def ConnectDefaultGenTrig(self,trigs=None,disconnect=False,partial=True):
    if not trigs:
      trigs = self.default_generation_triggers
    for trig in trigs:
      obj = self.getcont(trig)
      self.ConnectOneTrig(trig, obj, disconnect, partial)
      # another default generation trigger seems to be needed for the contentFlag
      # this is due to selecting a file from a previous job - the contentFlag appears to be set after the data object sends its dataChanged signal...
      if trig.startswith('F_SIGF') and not trig.endswith('_nonmtz'):
        self.ConnectOneTrig(trig+'.contentFlag', getattr(obj,'contentFlag'), disconnect, partial)

  @QtCore.Slot()
  def CellSpacegroup(self):
    ctrl, inp = self.container.controlParameters, self.container.inputData
    for var in ('CELL_A','CELL_B','CELL_C','CELL_D','CELL_E','CELL_F','SPACEGROUP'):
      self.getcont(var).setQualifier('allowUndefined',not inp.NON_MTZ)
      self.getWidget(var).validate()

  @QtCore.Slot(bool)
  def SADorMAD(self,init=False):
    ctrl, inp = self.container.controlParameters, self.container.inputData
    inp.F_SIGFanom.setQualifier(        'allowUndefined',bool(inp.NON_MTZ))
    inp.F_SIGFanom_nonmtz.setQualifier( 'allowUndefined',not inp.NON_MTZ)
    self.getWidget('F_SIGFanom').validate(), self.getWidget('F_SIGFanom_nonmtz').validate()
    inp.F_SIGFanom2.setQualifier(        'allowUndefined',not inp.MAD2 or bool(inp.NON_MTZ))
    inp.F_SIGFanom2_nonmtz.setQualifier( 'allowUndefined',not inp.MAD2 or not inp.NON_MTZ)
    inp.DNAME2.setQualifier(             'allowUndefined',not inp.MAD2 )
    self.getWidget('F_SIGFanom2').validate(),  self.getWidget('F_SIGFanom2_nonmtz').validate()
    self.getWidget('DNAME2').validate()
    inp.F_SIGFanom3.setQualifier(        'allowUndefined',not inp.MAD3 or bool(inp.NON_MTZ))
    inp.F_SIGFanom3_nonmtz.setQualifier( 'allowUndefined',not inp.MAD3 or not inp.NON_MTZ)
    inp.DNAME3.setQualifier(             'allowUndefined',not inp.MAD3 )
    self.getWidget('F_SIGFanom3').validate(),  self.getWidget('F_SIGFanom3_nonmtz').validate()
    self.getWidget('DNAME3').validate()
    inp.F_SIGFanom4.setQualifier(        'allowUndefined',not inp.MAD4 or bool(inp.NON_MTZ))
    inp.F_SIGFanom4_nonmtz.setQualifier( 'allowUndefined',not inp.MAD4 or not inp.NON_MTZ)
    inp.DNAME4.setQualifier(             'allowUndefined',not inp.MAD4 )
    self.getWidget('F_SIGFanom4').validate(),  self.getWidget('F_SIGFanom4_nonmtz').validate()
    self.getWidget('DNAME4').validate(),
    if not init:
      if inp.MAD2 or inp.MAD3 or inp.MAD4:
        #inp.EXPTYPE.set('MAD')
        # once the default pipeline is set in crank2, this can/should be moved there
        if self.container.controlParameters.USE_COMB:
          self.container.controlParameters.USE_COMB.set(False)
      else:
        self.SetSAD()
        self.InputNative()
    # this could be removed once it is possible to input "wildcard" programs to crank2
    if inp.MAD2 or inp.MAD3 or inp.MAD4 or (inp.EXPTYPE=='SIRAS' and inp.NATIVE):
      ctrl.PHAS_PROGRAM.set('bp3')
      ctrl.DMFULL_PHCOMB_PROGRAM.set('multicomb')
    else:
      ctrl.PHAS_PROGRAM.set('refmac')
      ctrl.DMFULL_PHCOMB_PROGRAM.set('refmac')

  @QtCore.Slot()
  def SetSAD(self):
    ctrl, inp = self.container.controlParameters, self.container.inputData
    #inp.EXPTYPE.set('SAD')
    if self.TASKNAME!='shelx':
      if inp.EXPTYPE=='SAD' and not ctrl.USE_COMB:
        ctrl.USE_COMB.set(True)
      if inp.EXPTYPE!='SAD' and ctrl.USE_COMB:
        ctrl.USE_COMB.set(False)

  @QtCore.Slot()
  def SSAD(self):
    ctrl, inp = self.container.controlParameters, self.container.inputData
    # this is in fact not only for S-SAD but also for all atoms that do not have wavelengths defined in crank2
    inp.DNAME.setQualifier('allowUndefined', inp.EXPTYPE!='MAD' and (self.TASKNAME=='shelx' or \
                             str(inp.ATOM_TYPE).strip().upper() not in ('SE','HG','AU','PT','ZN','BR')) )
    self.getWidget('DNAME').validate()

  @QtCore.Slot(bool)
  def InputNative(self,init=False):
    ctrl, inp = self.container.controlParameters, self.container.inputData
    inp.F_SIGFnative.setQualifier( 'allowUndefined',not inp.NATIVE or bool(inp.NON_MTZ) )
    inp.F_SIGFnative_nonmtz.setQualifier( 'allowUndefined',not inp.NATIVE or not inp.NON_MTZ)
    self.getWidget('F_SIGFnative').validate(),  self.getWidget('F_SIGFnative_nonmtz').validate()
    # once the default pipeline is set in crank2, this can/should be moved there
    #if inp.NATIVE and inp.ATOM_TYPE.isSet() and not init:
    #  # disabled as of now for Crank2 as both refmac and bp3 crash for SIRAS !  enabled again.lets see whether it was fixed...
    #  if inp.EXPTYPE=='SAD' and str(inp.ATOM_TYPE).strip()!='S'       : #and inp.SHELXCDE:
    #    inp.EXPTYPE.set('SIRAS')
    #    if self.container.controlParameters.USE_COMB:
    #      self.container.controlParameters.USE_COMB.set(False)
    #  if inp.EXPTYPE=='SIRAS' and str(inp.ATOM_TYPE).strip()=='S':
    #    self.SetSAD()
    #if not inp.NATIVE and inp.EXPTYPE=='SIRAS':
    self.SetSAD()

  def ClearSavedF(self,si):
    if not self.isEditable(): return
    inp = self.container.inputData
    self.ConnectDefaultGenTrig(disconnect=True)
    getattr(inp,'SAVED_FAVFILE'+si).unSet()
    if si!="_NATIVE":
      getattr(inp,'SAVED_FPMFILE'+si).unSet()
      if self.isEditable() and not self.job_started:
        getattr(inp,'DNAME'+si).unSet(),  getattr(inp,'WAVELENGTH'+si).unSet()
        getattr(inp,'FPRIME'+si).unSet(),  getattr(inp,'FDPRIME'+si).unSet()
        getattr(inp,'USER_WAVELENGTH'+si).set(False),  getattr(inp,'USER_DNAME'+si).set(False)
        getattr(inp,'USER_FPRIME'+si).set(False),  getattr(inp,'USER_FDPRIME'+si).set(False)
    self.ConnectDefaultGenTrig()

  def ToggleSADSIRASOption(self):
    return self.container.inputData.NATIVE and self.container.inputData.EXPTYPE in ('SAD','SIRAS')

  def ToggleHandDeterminationDoHand(self):
    return self.container.controlParameters.DO_HANDDET and self.basepipe.ToggleHandDetermination()

  def ToggleFree(self):
    return self.basepipe.CheckStartEnd('ref') or self.basepipe.CheckStartEnd('building')

  def ToggleNoUseComb(self):
    return self.basepipe.CheckStartEnd('building') and not self.basepipe.ToggleUseComb()
    # or self.container.inputData.INPUT_PARTIAL ) 

  def AutoSetI2DefaultPar(self, proc_parent, i2cont, par, par_recur=None):
    # tries to recursively find the corresponding default parameter in the proc_parent crank tree
    # the i2 par names (defined in the crank2 i2 xml) must match the crank2 names for this to work
    # returns 1 if parameter was set, 2 if it was not set and process not found, 0 if not set but process found
    if par_recur is None:
      par_recur=par
    proc_nick, sep, key = par_recur.lower().partition('_')
    # takes care of possible _ in the process nick
    if proc_nick not in [p.nick for p in proc_parent.processes]:
      proc_nick_add, sep, key = key.partition('_')
      proc_nick += '_'+proc_nick_add
    for proc in proc_parent.processes:
      if self.SI2P( proc, i2cont, par, key, proc_nick ):
        return 1
    else:
      for proc in proc_parent.processes:
        if proc.nick == proc_nick:
          if self.AutoSetI2DefaultPar(proc, i2cont, par, key)==1:
            return 1
          return 0
    return 2

  def SI2P(self, proc, i2cont, par, key, proc_nick):
    # only called internally by SetI2DefaultPar
    if proc.nick == proc_nick:
      if key == 'program':
        if proc.GetProg(supported=True):
          self.SetI2DefPar(i2cont, par, proc.GetProg(supported=True).nick)
          return 1
      elif key in proc.supported_params:
        self.SetI2DefPar(i2cont, par, proc.GetParam(key))
        return 1
    return 0

  def CheckUse(self, cont, param):
    # use this function for all parameters that can be set by default or by user
    if getattr(cont,param).isSet() and ( not self.has_cont_attr(cont,'USER_'+param) or getattr(cont,'USER_'+param) ):
      return True
    else:
      return False

  def SetI2DefPar(self,i2cont, par, value, allow_None=True, validate_widget=False, emulate=False):
    # return 1 if somethign was (re)set, 0 if nothing changed
    if not self.has_cont_attr(i2cont,'USER_'+par) or not getattr(i2cont,'USER_'+par) or not getattr(i2cont, par).isSet():
      if (value is not None or allow_None) and (str(value)!=str(getattr(i2cont, par))):
        con=getattr(i2cont, par)
        if not emulate:
          self.ConnectOneTrig(par,con,disconnect=True,partial=True,funct=self.setByUser)
          #if value==False:  value=0
          #elif value==True:  value=1
          if value=='unSet()':
            con.unSet()
          else:
            con.set(value)
          self.ConnectOneTrig(par,con,partial=True,funct=self.setByUser)
          if validate_widget and self.getWidget(par): 
            self.getWidget(par).validate()
        return 1
    return 0

  def CheckDataChanged(self,inp,lab):
    # checks whether data has changed and thus defaults need to be generated. returns 1 if yes, 0 if no
    if inp.NON_MTZ:
      return 1
    path, flag = str(getattr(inp,lab).fullPath), getattr(inp,lab).contentFlag
    # here we assume that input mtz will be converted to mini-mtz in ccp4 dir. without this we'd get
    # double generation and thus slow response. However if mini-mtz is not created then no generation
    # will take place - these two lines would have to be removed in such case!
    if not os.path.realpath(path).startswith(os.path.realpath(CCP4Utils.getTestTmpDir())) and \
       not os.path.realpath(path).startswith(os.path.realpath(PROJECTSMANAGER().getProjectDirectory(projectId=self.projectId()))):
      return 0
    try:
      flag = int(flag)
    except:
      flag = None
    # we don't need to generate defaults if there is no change in the path and flag
    if path==self.prev_data[lab][0] and flag==self.prev_data[lab][1]:
      return 0
    # saving the file path and content for the next check
    if path is not None and flag is not None:
      self.prev_data[lab]=(path,flag)
    return 1

  @QtCore.Slot(object)
  def setDefaultParameters(self,trig_by=None):
    from ....core import CCP4XtalData
    from ....core.CCP4PluginScript import CPluginScript
    ctrl = self.container.controlParameters
    inp  = self.container.inputData
    fsigf_supplied = True

    print('setDefaultParameters',trig_by)

    if trig_by is None:
      if self.isEditable():
        taskname=self.TASKNAME
        if inp.SHELX_SEPAR:
          self.TASKNAME='shelx'
        if self.TASKNAME=='shelx':
          inp.SHELX_SEPAR.set(True)
          self.ConnectDefaultGenTrig(trigs=['SHELXCDE','USE_COMB'],disconnect=True)
          inp.SHELXCDE.set(True)
          ctrl.FAEST_PROGRAM.set('shelxc')
          ctrl.SUBSTRDET_PROGRAM.set('shelxd')
          ctrl.USE_COMB.set(False)
          # different from the crank2 defaults - basically the (first) user is George here :)
          #ctrl.MBREF_BIGCYC.set(2)
          #ctrl.USER_MBREF_BIGCYC.set(True)
      self.SetPipeline(base_steps_only=not self.isEditable())
      if self.isEditable():
        if taskname.startswith('crank2_') and taskname.split('_')[1] in self.basepipe.base_steps:
          inp.START_PIPELINE.set( taskname.split('_')[1] )
        if taskname in ('crank2_comb_phdmmb','crank2_mbref'):
          inp.START_PIPELINE.set( 'building' )
        inp.F_SIGFnative_nonmtz.contentFlag.set( CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN )
        inp.F_SIGFanom_nonmtz.contentFlag.set( CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR )
        inp.F_SIGFanom2_nonmtz.contentFlag.set( CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR )
        inp.F_SIGFanom3_nonmtz.contentFlag.set( CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR )
        inp.F_SIGFanom4_nonmtz.contentFlag.set( CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR )
      if not ctrl.INITIAL_DEFAULTS_GEN:
        print('full setDefaultParameters not needed, quitting.')
        return
      else:
        ctrl.INITIAL_DEFAULTS_GEN.set(False)


    # we don't need to generate defaults if triggered by cell/spacegroup edit but not all supplied yet
    if trig_by and (trig_by.startswith('CELL_') or trig_by.startswith('SPACEGROUP')):
      if not inp.CELL_A or not inp.CELL_B or not inp.CELL_C or not inp.SPACEGROUP:
        return

    # clearing the SAVED F's info on new data input and checking whether the data change is relevant
    # (irrelevant changes may be triggered by i2 itself rather than user)
    if trig_by and trig_by.startswith('F_SIGFanom'):
      # we don't need to generate defaults if triggered by non_mtz objects but the non_mtz option is off (and the other way around)
      if ('_nonmtz' in trig_by and not inp.NON_MTZ) or ('_nonmtz' not in trig_by and inp.NON_MTZ):
        return
      si = trig_by[10]  if len(trig_by)>10 and trig_by[10].isdigit() else ''
      if not si or getattr(inp,'MAD'+si):
        if not self.CheckDataChanged(inp,'F_SIGFanom'+si):
          return
        self.ClearSavedF(si)
      else:
        return
    if trig_by and trig_by.startswith('F_SIGFnative'):
      if ('_nonmtz' in trig_by and not inp.NON_MTZ) or ('_nonmtz' not in trig_by and inp.NON_MTZ):
        return
      if getattr(inp,'NATIVE'):
        if not self.CheckDataChanged(inp,'F_SIGFnative'):
          return
        self.ClearSavedF('_NATIVE')
      else:
        return

    # check whether all the needed datasets were supplied
    for i in range(4):
      si = str(i+1)  if i  else ''
      sif = si+'_nonmtz'  if inp.NON_MTZ  else si
      if (not si or getattr(inp,'MAD'+si)) and (not getattr(inp,'F_SIGFanom'+sif).fullPath.isSet() or  \
                                                not getattr(inp,'F_SIGFanom'+sif).contentFlag.isSet()) :
        fsigf_supplied = False
        break
    if inp.NATIVE:
      if not inp.NON_MTZ and (not inp.F_SIGFnative.fullPath.isSet() or not inp.F_SIGFnative.contentFlag.isSet()):
        fsigf_supplied = False
      if inp.NON_MTZ and (not inp.F_SIGFnative_nonmtz.fullPath.isSet() or not inp.F_SIGFnative_nonmtz.contentFlag.isSet()):
        fsigf_supplied = False

    if self.isEditable() and not self.job_started and \
       ( (trig_by and trig_by.startswith('SEQIN')) or fsigf_supplied ) and \
       ( (trig_by and trig_by.startswith('F_SIGF')) or (inp.SEQIN.isSet() or not inp.INPUT_SEQUENCE) ):

      # check for multiple unassigned dnames and correct if it is the case
      nodn=0
      for dn,dnstr,use in [(inp.DNAME,'DNAME',True), (inp.DNAME2,'DNAME2',inp.MAD2), \
                       (inp.DNAME3,'DNAME3',inp.MAD3), (inp.DNAME4,'DNAME4',inp.MAD4)]:
        if use and (not dn.isSet() or not dn):
          if nodn:
            self.ConnectDefaultGenTrig(disconnect=True)
            self.SetI2DefPar( inp, dnstr, "unknown{}".format(nodn) )
            self.ConnectDefaultGenTrig()
          nodn+=1

      print("Starting generating defaults.")
      self.ConnectDefaultGenTrig(disconnect=True)
      try:
        workDir = PROJECTSMANAGER().jobDirectory(self.jobId(),subDir='TMP')
      except:
        workDir = None
      try:
        defaults = CPluginScript(dummy=True,workDirectory=workDir).makePluginObject(pluginName='crank2',reportToDatabase=False).process(self.container)
      except Exception as e:
        print("Exception caught when preparing the defaults: {}".format(e))
        print(traceback.print_exc())
        if 'MATTHEWS_COEF' in str(e) and ('The SYMM keyword' in str(e) or 'SYMMETRY OPERATOR ERROR' in str(e)):
          inp.SPACEGROUP.set(None)
      else:
        print("Defaults generation finished.",defaults)
        if defaults:
          try:
            # automatically setting defaults for the parameters named in the "crank2 way"
            for cont_str in ('controlParameters','inputData'):
              cont = getattr(self.container,cont_str)
              for par in cont._dataOrder:
                if not self.has_cont_attr(cont,'USER_'+par) or not getattr(cont,'USER_'+par):
                  if not self.AutoSetI2DefaultPar(defaults, cont, par) and par not in self.crank2_naming_exceptions:
                    print('WARNING : parameter {0} not found, typo suspected.'.format(par))

            # all the defaults for parameters that could not be set automatically should go here
            if defaults.GetParam('target'):
              self.SetI2DefPar(inp, 'EXPTYPE', defaults.GetParam('target'))
            first_proc = defaults.processes[0]  if defaults.processes[0].nick!='createfree'  else defaults.processes[1]
            if first_proc.inp.Get(has_solvent_content=True):
              self.SetI2DefPar(inp, 'SOLVENT_CONTENT', first_proc.inp.Get(has_solvent_content=True).solvent_content)
              self.SetI2DefPar(inp, 'MONOMERS_ASYM', first_proc.inp.Get(has_monomers_asym=True).monomers_asym)
            if first_proc.inp.Get(has_residues_mon=True):
              self.SetI2DefPar(inp, 'RESIDUES_MON', first_proc.inp.Get(has_residues_mon=True).residues_mon)
            if first_proc.inp.Get('model',typ=('substr','partial+substr'),has_num_atoms=True):
              #self.SetI2DefPar(inp, 'NUMBER_SUBSTRUCTURE', int(first_proc.inp.Get('model',typ=('substr','partial+substr'),has_num_atoms=True).exp_num_atoms)/int(inp.MONOMERS_ASYM))
              num_at=int(first_proc.inp.Get('model',typ=('substr','partial+substr'),has_num_atoms=True).exp_num_atoms)
              if num_at>0:
                self.SetI2DefPar(inp, 'NUMBER_SUBSTRUCTURE', num_at)
                inp.NUMBER_SUBSTRUCTURE.setQualifier('allowUndefined', True )
              # dealing with the obscure case of no S/SE atoms in the sequence (eg vapd)
              else:
                inp.NUMBER_SUBSTRUCTURE.setQualifier('allowUndefined', False )
              if self.getWidget('NUMBER_SUBSTRUCTURE'):
                self.getWidget('NUMBER_SUBSTRUCTURE').validate()
            #if defaults.GetProcess('substrdet') and defaults.GetProcess('substrdet').GetProg().IsKey('dsul'):
            #  self.SetI2DefPar(inp, 'SUBSTRDET_NUM_DSUL', int(defaults.GetProcess('substrdet').GetProg('shelxd').GetKey('dsul')) )
            if not ctrl.USE_COMB.isSet():
              ctrl.USE_COMB.set( True )
            if not ctrl.MB_PROGRAM.isSet() and self.basepipe.ToggleModelBuilding():
              if ctrl.USE_COMB:
                # it can happen that the program default is not assigned if eg cell exists but no spg
                if defaults.GetProcess('comb_phdmmb').GetProcess('mb') and defaults.GetProcess('comb_phdmmb').GetProcess('mb').GetProg(supported=True):
                  self.SetI2DefPar(ctrl, 'MB_PROGRAM', defaults.GetProcess('comb_phdmmb').GetProcess('mb').GetProg(supported=True).nick)
              elif defaults.GetProcess('mbref'): # ie assuming mbref
                if defaults.GetProcess('mbref').GetProcess('mb') and defaults.GetProcess('mbref').GetProcess('mb').GetProg(supported=True):
                  self.SetI2DefPar(ctrl, 'MB_PROGRAM', defaults.GetProcess('mbref').GetProcess('mb').GetProg(supported=True).nick)
            # shelx MIND second parameter needs to be "translated"
            if defaults.GetProcess('substrdet') and defaults.GetProcess('substrdet').IsParam('min_dist_symm_atoms'):
              self.SetI2DefPar(ctrl, 'SUBSTRDET_MIN_DIST_SYMM_ATOMS', True if defaults.GetProcess('substrdet').GetParam('min_dist_symm_atoms')<0.1 else False, validate_widget=True )
            # the first MIND parameter will show up positive but inputted as negative in the script
            if ctrl.SUBSTRDET_MIN_DIST_ATOMS.isSet() and ctrl.SUBSTRDET_MIN_DIST_ATOMS<0.0:
              self.SetI2DefPar(ctrl, 'SUBSTRDET_MIN_DIST_ATOMS', abs(ctrl.SUBSTRDET_MIN_DIST_ATOMS), validate_widget=True )

            def_inp = first_proc.inp
            mod = first_proc.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True)
            xn_nat = [o.GetCrystalName()  for o in def_inp.GetAll('fsigf',typ='average')  if 'native' in o.custom or o.GetCrystalName()=='native']
            xn = [o.GetCrystalName()  for o in def_inp.GetAll('fsigf',typ='plus')  if 'native' not in o.custom and o.GetCrystalName()!='native']
            # guess about substr atoms in native
            if inp.NATIVE and len(xn_nat)>0:
              self.SetI2DefPar(inp, 'SUBSTR_ATOMS_NATIVE', bool(defaults.inp.Get('model',typ=('substr','partial+substr'),xname=xn_nat[0])), validate_widget=True )
            #  self.SetI2DefPar(inp, 'SUBSTR_ATOMS_NATIVE', defaults.GetProcess('phdmmb').GetParam('substr_in_native'), validate_widget=True )
            dname_undef = mod.default_unknown if mod else 'Undefined'
            for i in range(4):
              si = str(i+1)  if i  else ''
              sif = si+'_nonmtz'  if inp.NON_MTZ  else si
              dn = getattr(inp,'DNAME'+si)  if getattr(inp,'DNAME'+si)  else dname_undef
              if dn in defaults.d_rename:
                dn = defaults.d_rename[dn]
              if getattr(inp,'F_SIGFanom'+sif).fullPath.isSet() and getattr(inp,'F_SIGFanom'+sif).contentFlag.isSet():
                if mod and inp.ATOM_TYPE.upper()==mod.GetAtomType() and not None in mod.Getfpfpp(mod.GetAtomType(),dn):
                  self.SetI2DefPar(inp, 'FPRIME'+si, '{0:.2f}'.format(mod.Getfpfpp(mod.GetAtomType(),dn)[0]) )
                  self.SetI2DefPar(inp, 'FDPRIME'+si, '{0:.2f}'.format(mod.Getfpfpp(mod.GetAtomType(),dn)[1]) )
                if first_proc.inp.Get('fsigf',dname=dn,has_wavel=True):
                  self.SetI2DefPar(inp, 'WAVELENGTH'+si, first_proc.inp.Get('fsigf',dname=dn,has_wavel=True).wavel)
                # assigning the F's to speed up the future defaults generation.
                if def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='plus',col='f'):
                  getattr(inp,'SAVED_FPMFILE'+si).set(  def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='plus',col='f').GetFileName() )
                  getattr(inp,'SAVED_FAVFILE'+si).set(  def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='average',col='f').GetFileName() )
                  getattr(inp,'SAVED_FPLUS'+si).set(    def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='plus',col='f').GetLabel('f') )
                  getattr(inp,'SAVED_SIGFPLUS'+si).set( def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='plus',col='f').GetLabel('sigf') )
                  getattr(inp,'SAVED_FMIN'+si).set(     def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='minus',col='f').GetLabel('f') )
                  getattr(inp,'SAVED_SIGFMIN'+si).set(  def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='minus',col='f').GetLabel('sigf') )
                  getattr(inp,'SAVED_FAVER'+si).set(    def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='average',col='f').GetLabel('f') )
                  getattr(inp,'SAVED_SIGFAVER'+si).set( def_inp.Get('fsigf',xname=xn[0],dname=dn,typ='average',col='f').GetLabel('sigf') )
                if not si and xn_nat and def_inp.Get('fsigf',xname=xn_nat[0],typ='average',col='f'):
                  getattr(inp,'SAVED_FAVFILE_NATIVE').set( def_inp.Get('fsigf',xname=xn_nat[0],typ='average',col='f').GetFileName() )
                  getattr(inp,'SAVED_FAVER_NATIVE').set(   def_inp.Get('fsigf',xname=xn_nat[0],typ='average',col='f').GetLabel('f') )
                  getattr(inp,'SAVED_SIGFAVER_NATIVE').set(def_inp.Get('fsigf',xname=xn_nat[0],typ='average',col='f').GetLabel('sigf') )
                # reassign DNAME 
                if mod and mod.Getfpfpp(mod.GetAtomType(),dn,return_guessed_dtype=True)[4]:
                  self.SetI2DefPar( inp,'DNAME'+si, mod.Getfpfpp(mod.GetAtomType(),dn,return_guessed_dtype=True)[4], validate_widget=True )
            # assign cell & spacegroup
            if inp.NON_MTZ:
              asc,ass=False,False
              for fsigf in def_inp.GetAll('fsigf',xname=xn[0]) if xn else [] + def_inp.GetAll('fsigf',xname=xn_nat[0]) if xn_nat else []:
                if fsigf.cell:
                  for i,c in enumerate(('CELL_A','CELL_B','CELL_C','CELL_D','CELL_E','CELL_F')):
                    if not self.SetI2DefPar(inp,c,fsigf.cell[i], emulate=True):
                      break
                  else:
                    for i,c in enumerate(('CELL_A','CELL_B','CELL_C','CELL_D','CELL_E','CELL_F')):
                      self.SetI2DefPar(inp,c,format(float(fsigf.cell[i]),'.1f'), validate_widget=True)
                  asc=True
                if fsigf.spgr:
                  self.SetI2DefPar(inp,'SPACEGROUP',fsigf.spgr, validate_widget=True)
                  ass=True
                if asc and ass:
                  break

          except Exception as e:
            print("Exception caught when setting the defaults: {}".format(e))
            print(traceback.print_exc())

          print("Defaults assigned.")

      self.ConnectDefaultGenTrig()

  @QtCore.Slot(str)
  def setByUser(self, param=None):
    if self.getcont('USER_'+param) is not None:
      print('Parameter {0} set by user'.format(param))
      self.getcont('USER_'+param).set( True )

  def taskValidity(self):
    rv = CCP4ErrorHandling.CErrorReport()
    ctrl, inp = self.container.controlParameters, self.container.inputData
    if not inp.NON_MTZ:
      for i in range(4):
        si = str(i+1)  if i  else ''
        if (not si or getattr(inp,'MAD'+si)) and (not getattr(inp,'F_SIGFanom'+si).contentFlag.isSet() or \
                                                  getattr(inp,'F_SIGFanom'+si).contentFlag not in (0,1,2)):
          rv.append(self.__class__,' Wrong input data F_SIGFanom\n',details='Input file {}   does not contain anomalous (I+,I- or F+,F-) data.'.format(getattr(inp,'F_SIGFanom'+si).fullPath),stack=False)
          rv[-1]['description']="No anomalous data inputted for F_SIGFanom{}.".format(si)
    if inp.SUBSTRDET_NUM_DSUL.isSet() and inp.NUMBER_SUBSTRUCTURE and inp.SUBSTRDET_NUM_DSUL*2>inp.NUMBER_SUBSTRUCTURE:
      rv.append(self.__class__,' Too many disulphides inputted\n',details='{} sulpher atoms may form at most {} disulphides, however, {} disulphides inputted.'.format(
          inp.NUMBER_SUBSTRUCTURE,int(inp.NUMBER_SUBSTRUCTURE)/2,inp.SUBSTRDET_NUM_DSUL), stack=False)
      rv[-1]['description']="Number of disulphides larger than number of sulphers divided by two."
    if self.TASKNAME!='shelx' and (not inp.FPRIME or not inp.FDPRIME):
      rv.append(self.__class__,' f\' or f\'\' not known\n',details='Anomalous scattering coefficients f\', f\'\' should be directly inputted especially if they are known from a fluorescence scan. Otherwise data type or wavelength input can be used to calculate their theoretical value.', stack=False)
      rv[-1]['description']="Anomalous scattering coefficients are needed by Crank2."
    if len(rv)>0:
      self.job_started=False
    return rv

  # defining job_started so that DNAME etc is not cleared on i2 internal data file path changes etc
  def fix(self):
    self.job_started=True
    return CCP4ErrorHandling.CErrorReport()

def PrepCont(cont,name='crank2'):
  # helper function preparing a new container from an existing one
  from ....core import CCP4Container
  from ....core import CCP4File
  cont_new = CCP4Container.CContainer(name=name)
  cont_new.__dict__['header'] = CCP4File.CI2XmlHeader()
  cont_new.__dict__['header'].pluginName = 'crank2' #name
  for ic in ('inputData','controlParameters'):
    cont_new.addObject( getattr(cont,ic), reparent=True )
  cont_new.controlParameters.INITIAL_DEFAULTS_GEN.set(True)
  #cont_new.inputData.END_PIPELINE.set('ref')
  cont_new.inputData.END_PIPELINE.set('building')
  cont_new.inputData.SHELX_SEPAR.set(name=='shelx'), cont_new.inputData.SHELXCDE.set(name=='shelx')
  return cont_new

def whatNext(jobId,childTaskName,childJobNumber,projectName):
  whatnext = []
  cont = PROJECTSMANAGER().getJobParams(jobId)
  if cont.outputData.XYZOUT.isSet():
    whatnext.extend(['coot_rebuild','prosmart_refmac'])
  if cont.outputData.FPHOUT_HL.isSet():
    whatnext.append(['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml'])
    whatnext.append('arp_warp_classic')
    if not cont.outputData.XYZOUT.isSet():
      whatnext.extend(['parrot'])
  if str(cont.inputData.EXPTYPE)!='MAD' or not cont.outputData.XYZOUT.isSet():
    jobdir = PROJECTSMANAGER().jobDirectory(jobId)
    xml_new, xml_sh = os.path.join(jobdir,'i2_crank2_next.xml'), os.path.join(jobdir,'i2_shelx_next.xml')
    if cont.outputData.XYZOUT_SUB_RES.isSet() and cont.outputData.XYZOUT_SUBSTR.isSet():
      cont_sh=PrepCont(cont,'shelx')
      cont_sh.inputData.XYZIN_SUB_RES.set( cont.outputData.XYZOUT_SUB_RES )
      cont_sh.inputData.XYZIN_SUB.set( cont.outputData.XYZOUT_SUBSTR )
      cont_sh.inputData.START_PIPELINE.set('phdmmb')
      cont_sh.saveDataToXml(fileName=xml_sh)
      whatnext.append(['shelx',xml_sh])
    cont_new=PrepCont( PROJECTSMANAGER().getJobParams(jobId) )
    cont_new.controlParameters.USE_COMB.set(True)
    if cont.outputData.F_SIGFanom_OUT.isSet():
      cont_new.inputData.F_SIGFanom.set( cont.outputData.F_SIGFanom_OUT )
      if bool(cont_new.inputData.NON_MTZ):
        cont_new.inputData.NON_MTZ.set(False)
    if cont.outputData.XYZOUT_SUBSTR.isSet():
      cont_new.inputData.XYZIN_SUB.set( cont.outputData.XYZOUT_SUBSTR )
      cont_new.inputData.START_PIPELINE.set('refatompick')
    if cont.outputData.FPHOUT_HL.isSet():
      cont_new.inputData.INPUT_PHASES.set(True)
      cont_new.inputData.FPHIN_HL.set( cont.outputData.FPHOUT_HL )
      cont_new.inputData.START_PIPELINE.set('dmfull')
    if cont.outputData.XYZOUT.isSet():
      cont_new.inputData.INPUT_PARTIAL.set(True)
      cont_new.inputData.XYZIN.set( cont.outputData.XYZOUT )
      cont_new.inputData.START_PIPELINE.set('building')
    if str(cont.inputData.EXPTYPE)!='MAD' or not cont.outputData.XYZOUT.isSet():
      cont_new.saveDataToXml(fileName=xml_new)
    whatnext.append(['crank2',xml_new])
    #whatnext.append([cont_new.__dict__['header'].pluginName,xml_new])
  return whatnext


def exportMtzColumnLabels(jobId=None,paramNameList=None,sourceInfoList=[]):
  colLabels = { 'FPHOUT_HL':'FPHOUT_HL', 'FPHOUT_DIFF':'FPHOUT_DIFF', 'FPHOUT_2FOFC':'FPHOUT_2FOFC', 'FPHOUT_DIFFANOM':'FPHOUT_DIFFANOM' }
  #print 'CTaskCrank2.exportMtzColumnLabel',jobId,paramNameList,sourceInfoList
  from ....core import CCP4Container
  paramsFile = PROJECTSMANAGER().makeFileName(jobId = jobId,mode='PARAMS')
  #print 'CTaskCrank2.exportMtzColumnLabel paramsFile',paramsFile
  c = CCP4Container.CContainer()
  c.loadDataFromXml(paramsFile)
  ret = []
  indx = 0
  for  paramName in paramNameList:
    if paramName.startswith('F_SIGF'):
      if paramName.startswith('F_SIGFanom'):
        dParam = 'DNAME' + paramName[-1] if paramName[-1] in ['2','3','4'] else 'DNAME'
        label = c.inputData.get(dParam).__str__() if c.inputData.get(dParam) else 'anom'
      else:
        label = 'native'
      ret.append('crank_'+label)
    elif colLabels.get(paramName,None) is not None:
      ret.append('crank_'+colLabels[paramName])
    else:
      indx += 1
      ret.append('crank_'+str(indx))
  return ret
