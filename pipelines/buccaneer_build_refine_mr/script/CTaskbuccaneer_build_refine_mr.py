"""
     tasks/buccaneer_build_refine/CTaskbuccaneer_build_refine.py: CCP4 GUI Project
     Copyright (C) 2011 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
     Liz Potterton Dec 2011 - copy and revise buccaneer gui
     Jon Agirre 2014 - Add new MR options, revise the gui and simplify it
"""

from baselayer import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets


def whatNext(jobId,childTaskName,childJobNumber,projectName):
  #print 'buccaneer_build_refine_mr.whatNext',projectName
  if childJobNumber is None:
    return ['coot_rebuild','edstats','prosmart_refmac','modelcraft']
  elif childTaskName == 'refmac':
    return  ['coot_rebuild','edstats','prosmart_refmac','modelcraft']
  else:
    return [['prosmart_refmac','Refine this model before continuing']]


# Subclass CTaskWidget to give specific task window
class CTaskbuccaneer_build_refine_mr(CCP4TaskWidget.CTaskWidget):
  TASKNAME = 'buccaneer_build_refine_mr'
  TASKVERSION = 0.0
  TASKMODULE='deprecated'
  TASKTITLE='Autobuild protein - BUCCANEER'
  SHORTTASKTITLE='BUCCANEER'
  TASKLABEL = 'bucanr'
  DESCRIPTION = 'Iterations of model building (Buccaneer) and refinement (Refmac5, Prosmart and Coot)'
  MGDISPLAYFILES = ['XYZOUT']
  WHATNEXT = ['coot_rebuild','prosmart_refmac','modelcraft', 'edstats']
  PROGRAMHELP = ['cbuccaneer','refmac5']
  RANK=1
  EXPORTMTZPARAMS = [ 'F_SIGF','FPHIOUT', ['DIFFPHIOUT' ,'diff' ] , 'ABCDOUT' ]

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)


  def drawContents(self):

    self.setProgramHelpFile('buccaneer')

    # Trying to handle if task is cloned and ABCD already set - then need to set ABCD_MR or ABCD_EXP
    if self.container.inputData.ABCD.isSet():
      if self.container.controlParameters.BUCCANEER_PHSIN_TYPE == 'mr' :
        self.container.inputData.ABCD_MR.set(self.container.inputData.ABCD)
      else:
        self.container.inputData.ABCD_EXP.set(self.container.inputData.ABCD)
      self.container.inputData.ABCD.unSet()
    self.container.inputData.ABCD_EXP.setQualifier ( 'allowUndefined', False )
    

# the input data tab starts here    
    
    folder = self.openFolder(folderFunction='inputData',title='Input Data')
    self.createLine( [ 'label', 'Build model with phases coming from  ', 'widget', '-guiMode', 'radio', 'BUCCANEER_PHSIN_TYPE' ] )

    self.createLine( [ 'subtitle', 'Input the molecular replacement model used to phase the data', 'A ' +\
                       '<i>refined</i> molecular replacement solution provides best initial phases for the ' +\
                       'autobuilding procedure; jelly body refinement will be done if the supplied model ' +\
                       'has not been refined previously.' ], toggle=['BUCCANEER_PHSIN_TYPE', 'open', 'mr' ] )
    self.openSubFrame(frame=[True], toggle= [ 'BUCCANEER_PHSIN_TYPE', 'open', 'mr' ] )
    self.createLine( [ 'widget','BUCCANEER_MR_MODE_XYZIN' ] )
    self.createLine( [ 'label', 'This model will be used to place and name chains, and','widget','BUCCANEER_MR_MODE' ] )
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Select experimental data', 'Looking for a way to input your pre-computed phases? You can find it under <i>Advanced Buccaneer Options</i>.' ]
                     , toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'mr'] )
    self.createLine( [ 'subtitle', 'Select experimental data', 'Experimental phases will no longer be used as restraints if R<sub>work</sub> &lt; ' + str( self.getContainer().controlParameters.BUCCANEER_MLHL_RWORK_LIMIT ) ]
                     , toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'experimental'] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'widget', 'F_SIGF' ] )
    self.createLine( [ 'widget', 'ABCD_EXP' ], toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'experimental'] )
    
    self.createLine( [ 'widget', 'FREERFLAG' ] )
    self.closeSubFrame ( )

    self.createLine( [ 'subtitle', 'Enter the AU content containing the structure sequence(s)' ] )
    self.openSubFrame ( frame=[True] )
    
    self.createLine( [ 'widget','ASUIN' ] )
    self.closeSubFrame ( )

    self.createLine( [ 'tip', 'Your model should contain confidently-built metals, ligands and/or protein parts', 'widget', 'XYZIN_MODE', 'label', 'Start from a partially built model' ] )
    self.container.controlParameters.XYZIN_MODE.dataChanged.connect( self.updateRequirements )
    self.container.controlParameters.BUCCANEER_PHSIN_TYPE.dataChanged.connect( self.updateRequirements )
    
    self.openSubFrame( frame=[True], toggle=[ 'XYZIN_MODE', 'open', [True] ] )
    self.createLine( [ 'widget','XYZIN' ] )
    self.createLine( [ 'tip', 'Syntax reminder: /ch/res/at/:radius e.g. /A/*/*:2.0 for all atoms in chain A. You need not specify the radius.', 'label', 'Optionally keep part of the known structure fixed ','widget', 'KNOWN_STRUCTURE' ] )
    self.closeSubFrame()
    self.closeFolder()

# a control parameter tab with basic options starts here

    folder = self.openFolder( folderFunction='controlParameters', title='Options' , drawFolder= self.drawOptions )

# a control parameter tab with more advanced options starts here

    folder = self.openFolder(folderFunction='controlParameters',title='Advanced Buccaneer Options',drawFolder=self.drawAdvancedBuccaneer )
    

# a control parameter tab with coot and refmac options starts here     

    folder = self.openFolder(folderFunction='controlParameters',title='Refinement options',drawFolder=self.drawRefinement)


# the last tab, with the rarely used options for choosing reference structures

    folder = self.openFolder(folderFunction='controlParameters',title='Reference structures',drawFolder=self.drawReference)
    

    #for w in self.findChildren(CCP4TaskWidget.CStackedWidget):
    #  print '   ',w,w.controlVar

    self.updateRequirements()
    self.updateViewFromModel()
    

  @QtCore.Slot()
  def updateRequirements ( self ) :
    
    if self.container.controlParameters.BUCCANEER_PHSIN_TYPE == 'mr' :
      self.container.inputData.BUCCANEER_MR_MODE_XYZIN.setQualifier ( 'allowUndefined', False )
      self.container.inputData.BUCCANEER_MR_MODE_XYZIN.setQualifier ( 'toolTip', 'A molecular replacement model is required' )
    else :
      self.container.inputData.BUCCANEER_MR_MODE_XYZIN.setQualifier ( 'allowUndefined', True )
    self.getWidget ( 'BUCCANEER_MR_MODE_XYZIN' ).validate ( )
    
    if self.container.controlParameters.XYZIN_MODE :
      self.container.inputData.XYZIN.setQualifier ( 'allowUndefined', False )
      self.container.inputData.XYZIN.setQualifier ( 'toolTip', 'An initial atomic model is required' )
    else :
      self.container.inputData.XYZIN.setQualifier ( 'allowUndefined', True )
    self.getWidget ( 'XYZIN' ).validate ( )
    
    return

  def fix(self):
    self.container.inputData.ABCD.unSet()
    if str(self.container.controlParameters.BUCCANEER_PHSIN_TYPE) == 'mr' :
      self.container.inputData.ABCD.set(self.container.inputData.ABCD_MR)
    else:
      self.container.inputData.ABCD.set(self.container.inputData.ABCD_EXP)
    self.container.inputData.ABCD_EXP.setQualifier ( 'allowUndefined', True )
    self.getWidget('ABCD_EXP').validate()
    self.container.inputData.ABCD_EXP.unSet()
    self.container.inputData.ABCD_MR.unSet()
    return CCP4TaskWidget.CTaskWidget.fix(self)


  def drawOptions(self):
    self.createLine ( [ 'subtitle', 'Pipeline control', '<b>Warning</b>: default values work best for most situations.' ] )
    self.openSubFrame ( frame = [True] )
    self.createLine( [ 'widget', 'ITERATIONS', 'label', 'Pipeline iterations'] )
    self.createLine( [ 'widget', 'STOP_AUTOMATICALLY', 'label', 'Stop automatically if there is no improvement'] )
    self.createLine( [ 'widget', 'BUCCANEER_CYCLES', 'label', 'initial buccaneer cycles' ] )
    self.createLine( [ 'label', 'For each cycle perform:' ] )
    self.createLine( [ 'widget', 'BUCCANEER_CYCLES_NEXT', 'label', 'buccaneer cycles'] )
    self.createLine( [ 'widget', 'REFMAC_CYCLES', 'label', 'refinement cycles'] )
    self.createLine( [ 'label', 'Extra steps:'])
    self.createLine( [ 'widget', 'CHAIN_PRUNE', 'label', 'Delete badly scoring chains (&#8804;20 residues) at the end of each cycle.'])
    self.createLine( [ 'widget', 'FULL_PRUNE', 'label', 'Delete badly scoring chains, residues and side chains at the start of each cycle. This usually helps at high resolution (&lt;2.5A).'])
    self.closeSubFrame ( )

    self.createLine ( [ 'subtitle', 'Buccaneer options', '<b>Warning</b>: default values work best for most situations.' ] )
    self.openSubFrame ( frame = [True] )
    self.createLine( [ 'widget', 'BUCCANEER_ANISOTROPY_CORRECTION','label','Apply anisotropy correction' ] )
    self.createLine( [ 'widget', 'BUCCANEER_BUILD_SEMET', 'label','Build methionines as selenomethionine' ] )
    self.createLine( [ 'advice', 'The three-letter code MSE will be used instead of MET' ], toggle=[ 'BUCCANEER_BUILD_SEMET', 'open', [True] ] )
    
    self.createLine( [ 'widget', 'BUCCANEER_FAST','label','Use fast mode - recommended' ] )
    self.closeSubFrame ( )

  def drawAdvancedBuccaneer(self ):
    
    self.createLine( [ 'subtitle', 'Filtering', 'Buccaneer calculates the mean B<sub>factor</sub> distribution for all residues, and leaves out those showing values higher than a chosen &sigma; level' ] )
    self.openSubFrame ( frame=[True])

    self.createLine( [ 'label', 'Filter the molecular replacement model by a', 'widget', 'BUCCANEER_MR_MODE_SIGMA', 'label', '&sigma; threshold'], toggle=['BUCCANEER_PHSIN_TYPE', 'open', 'mr' ] )
    self.createLine( [ 'label', 'Filter the model between cycles by a', 'widget', 'BUCCANEER_MODEL_SIGMA', 'label', '&sigma; threshold' ] )
    self.createLine( [ 'tooltip', '(Using data beyond 2A is slow and seldom helps)', 'label','Leave out reflection data beyond a','widget','BUCCANEER_RESOLUTION', 'label', '&#8491; resolution limit'] )
    self.createLine( [ 'widget' , 'BUCCANEER_USE_FREER', 'label', 'exclude Free R reflections from map calculation' ] )
    self.closeSubFrame ( )

    self.createLine( [ 'subtitle', 'Input a heavy atom model to aid with sequencing', 'You may input, for instance, Selenium atoms to aid in placing SeMet residues.' ] )
    self.openSubFrame ( frame = [True] )
    self.createLine( [ 'widget', 'XYZIN_SEQ' ] )
    self.closeSubFrame ( )
    
    self.createLine( [ 'subtitle', 'Pre-computed phases', 'If you have already refined your model, you may input your phase probabilities here.' ], toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'mr'] )
    self.openSubFrame ( frame=[True], toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'mr'])
    self.createLine( [ 'widget', 'ABCD_MR' ] )
    self.closeSubFrame ( ) 
    self.createLine( [ 'widget','BUCCANEER_FIX_POSITION','label','Build the new model in the same place as the input model' ] )
    self.createLine( [ 'label','Residue name for unsequenced residues','widget','BUCCANEER_NEW_RESIDUE_NAME' ] )
  
    import multiprocessing, os
    self.getContainer().controlParameters.BUCCANEER_JOBS = int ( os.getenv ( 'OMP_NUM_THREADS', multiprocessing.cpu_count() ) )
    self.createLine( [ 'tooltip',  'The default value shows the detected number of cores', 'label', 'Use up to', 'widget', 'BUCCANEER_JOBS', 'label', 'CPU cores whenever possible' ] ) 
    self.createLine( [ 'widget','BUCCANEER_CLEANUP','label','Clean up intermediate map files as job progresses' ] )


  def drawRefinement(self):
    
    self.createLine( [ 'subtitle', 'Real space refinement', 'A number of operations are available to be performed between the building and reciprocal-space refinement steps' ] )
    
    self.openSubFrame( frame=[True] )
    self.createLine( [ 'label', 'Perform ','widget','COOT_REALSPACE_OPERATION' , 'label', ' on the built model'] )   
    self.createLine( [ 'label', '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;when R<sub>work</sub> is ', 'widget', 'BUCCANEER_RSR_RWORK_LIMIT', 'label', ' or lower' ], toggle = [ 'COOT_REALSPACE_OPERATION','open', [ 'coot_stepped_refine','coot_fit_residues', 'coot_script_lines', 'coot_add_waters' ] ]  )
    self.createLine( [ 'widget', 'USERAMA', 'label', 'Use ramachandran restraints' ], toggle = [ 'COOT_REALSPACE_OPERATION','open', [ 'coot_stepped_refine','coot_fit_residues' ] ]  )
    
    self.createLine( [ 'label','Script to run' ], toggle = [ 'COOT_REALSPACE_OPERATION','open', [ 'coot_script_lines' ] ] )
    self.createLine( [ 'widget', '-guiMode','multiLine','SCRIPT' ], toggle = [ 'COOT_REALSPACE_OPERATION','open', [ 'coot_script_lines' ] ] )
    
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Reciprocal space refinement', '<b>Note:</b></br>Default options need to be changed very rarely. Use with caution' ] )
    
    self.openSubFrame( frame=[True] )
    self.createLine( [ 'tip', 'This option will be turned off automatically if the generated model is very fragmented', 'widget', 'REFMAC_LOCAL_NCS', 'label', 'Exploit local non-crystallographic symmetry' ] )
    #self.createLine( [ 'advice', 'This option will be turned off automatically if the generated model is very fragmented' ], toggle=['REFMAC_LOCAL_NCS', 'open', [True] ] )
    self.createLine( [ 'widget', 'REFMAC_EXP_USEPHI', 'label', 'Use experimental phases as restraints in refinement' ], toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'experimental' ] )
    self.createLine( [ 'label', '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;but drop them if R<sub>work</sub> goes below', 'widget', 'BUCCANEER_MLHL_RWORK_LIMIT' ], toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'experimental' ] )
    self.createLine( [ 'widget', 'REFMAC_MR_USEPHI', 'label', 'Use molecular replacement phases as restraints in refinement' ], toggle=[ 'BUCCANEER_PHSIN_TYPE', 'open', 'mr' ] )
    self.createLine( [ 'widget', 'USE_PROSMART', 'label', 'Use a higher-resolution restraint target' ] )
    self.createLine( [ 'widget', 'TARGET' ], toggle=[ 'USE_PROSMART', 'open', [True]] )
    self.createLine( [ 'widget', 'DICT'   ] )
    self.createLine( [ 'widget', '-guiMode','multiLine','EXTRAREFMACKEYWORDS' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Other', 'Additional MR model refinement' ] )
    self.openSubFrame( frame=[True] )
    self.createLine( [ 'widget', 'USE_SHIFTFIELD','label','Use shift field refinement of an initial MR model' ] )
    self.closeSubFrame()

  def drawReference(self):
    self.createLine(  [ 'subtitle', 'Reference structures', 'You should normally let Buccaneer to choose reference structures'] )
    self.openSubFrame ( frame=[True] )
    self.createLine(  [ 'widget', 'F_SIGF_REF' ] )
    self.createLine(  [ 'widget', 'ABCD_REF' ] )
    self.createLine(  [ 'widget', 'XYZIN_REF' ] )
    self.closeSubFrame ()
