"""
    servalcat_pipe_gui.py: CCP4 GUI Project
    Copyright (C) 2024 University of Southampton, MRC LMB Cambridge

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

from PySide2 import QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from core import CCP4XtalData
from pipelines.import_merged.script.dybuttons import ChoiceButtons
import os
import shutil
import gemmi


def whatNext(jobId=None,childTaskName=None,childJobNumber=None,projectName=None):
    import os
    from core import CCP4Modules, CCP4Utils, CCP4File, CCP4Container, CCP4Data, CCP4PluginScript
    jobStatus = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId,'status')
    if jobStatus == 'Unsatisfactory':
        returnList = ['LidiaAcedrgNew', 'servalcat_pipe']
    else:
        returnList = ['servalcat_pipe', 'coot_rebuild', 'modelcraft']
    return returnList

class Cservalcat_pipe(CCP4TaskWidget.CTaskWidget):

  TASKTITLE='Refinement - Servalcat'
  SHORTTASKTITLE='Servalcat'
  DESCRIPTION='Refinement against diffraction data or cryo-EM SPA map with optional restraints from ProSmart and/or MetalCoord'
  TASKNAME = 'servalcat_pipe'
  TASKLABEL = 'servalcat_pipe'
  TASKVERSION = 0.0
  TASKMODULE = 'refinement'
  MGDISPLAYFILES = ['XYZOUT','FPHIOUT','DIFFPHIOUT']
  AUTOPOPULATEINPUT = True
  RANK=1

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def ToggleWeightAdjustRmszAvailable(self):
    return str(self.container.controlParameters.WEIGHT_OPT) == 'AUTO' and \
      not self.container.controlParameters.WEIGHT_NO_ADJUST

  def ToggleRigidModeOn(self):
    return str(self.container.controlParameters.REFINEMENT_MODE) == 'RIGID'

  def ToggleRigidModeOff(self):
    return not str(self.container.controlParameters.REFINEMENT_MODE) == 'RIGID'

  def ToggleRestraintsOn(self):
    return (not bool(self.container.controlParameters.UNRESTRAINED)) and (not bool(self.container.controlParameters.FIX_XYZ)) and (not bool(self.container.controlParameters.JELLY_ONLY))

  def ToggleRestraintsOff(self):
    return bool(self.container.controlParameters.UNRESTRAINED) or bool(self.container.controlParameters.FIX_XYZ) or bool(self.container.controlParameters.JELLY_ONLY)

  def ToggleMetalCoordGenerate(self):
    return bool(self.container.metalCoordPipeline.RUN_METALCOORD) and str(self.container.metalCoordPipeline.GENERATE_OR_USE) == "GENERATE"

  def ToggleMetalCoordGenerateAdvanced(self):
    return bool(self.container.metalCoordPipeline.RUN_METALCOORD) and str(self.container.metalCoordPipeline.GENERATE_OR_USE) == "GENERATE" and bool(self.container.metalCoordPipeline.TOGGLE_ADVANCED)

  def ToggleMetalCoordUse(self):
    return bool(self.container.metalCoordPipeline.RUN_METALCOORD) and str(self.container.metalCoordPipeline.GENERATE_OR_USE) == "USE"

  def ToggleJellyOn(self):
    return (not bool(self.container.controlParameters.UNRESTRAINED)) and (not bool(self.container.controlParameters.FIX_XYZ))

  def ToggleJellyOff(self):
    return bool(self.container.controlParameters.UNRESTRAINED) or bool(self.container.controlParameters.FIX_XYZ)

  def ToggleFSIGFOrISIGIavailable(self):
     return (bool(self.container.controlParameters.HKLIN_IS_I_SIGI) and str(self.container.controlParameters.MERGED_OR_UNMERGED) == "merged")

  def ToggleTwinSuboptimal(self):
    if str(self.container.controlParameters.MERGED_OR_UNMERGED) == "unmerged":
       return False
    return (not self.container.controlParameters.HKLIN_IS_I_SIGI or self.container.controlParameters.F_SIGF_OR_I_SIGI == "F_SIGF")

  def ToggleOccComplete(self):
    if self.container.controlParameters.OCCUPANCY_GROUPS and self.container.controlParameters.OCCUPANCY_COMPLETE:
        return True
    return False

  def ToggleOccIncomplete(self):
    if self.container.controlParameters.OCCUPANCY_GROUPS and self.container.controlParameters.OCCUPANCY_INCOMPLETE:
        return True
    return False

  def drawContents(self):
    self.setProgramHelpFile('servalcat')
    if self.container.metalCoordPipeline.LIGAND_CODES_AVAILABLE:
        self.monomersWithMetals = self.container.metalCoordPipeline.LIGAND_CODES_AVAILABLE
        self.container.metalCoordPipeline.LIGAND_CODES_SELECTED = self.container.metalCoordPipeline.LIGAND_CODES_AVAILABLE
    else:
        self.monomersWithMetals = []
    #-  --------------------          --------------------          --------------------
    self.openFolder(folderFunction='inputData',title='Input Data')
    self.hklinChanged()

    self.createLine( [ 'subtitle', 'Main inputs' ])
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', '-browseDb', True, 'XYZIN' ])
    self.closeSubFrame()
    if self.isEditable():
       self.container.inputData.XYZIN.dataChanged.connect( self.modelChanged)
    # self.createLine( [ 'label', 'Experimental data type:', 'widget', 'DATA_METHOD' ])

    self.openSubFrame(toggle = ['DATA_METHOD', 'open', [ 'xtal' ] ] )
    self.createLine( [ 'widget', '-browseDb', True, 'HKLIN' ], toggle = ['MERGED_OR_UNMERGED', 'open', [ 'merged' ] ] )
    self.createLine( [ 'widget', '-browseDb', True, 'HKLIN_UNMERGED' ], toggle = ['MERGED_OR_UNMERGED', 'open', [ 'unmerged' ] ] )
    self.createLine( [ 'label', '<i>Warning: Unmerged data should be scaled.</i>' ], toggle = ['MERGED_OR_UNMERGED', 'open', [ 'unmerged' ] ] )
    self.container.inputData.HKLIN.dataChanged.connect( self.hklinChanged )
    self.createLine( [ 'label', 'Diffraction data are', 'widget', 'MERGED_OR_UNMERGED'] )
    self.createLine( [ 'label', 'Refinement against <b>amplitudes</b>.'], toggle = ['HKLIN_IS_I_SIGI', 'open', [ False ] ] )
    self.createLine( [ 'label', 'Refinement against <b>intensities</b>.'], toggle = ['MERGED_OR_UNMERGED', 'open', [ 'unmerged' ] ] )
    self.createLine( [ 'label', 'Refinement against', 'widget', 'F_SIGF_OR_I_SIGI'], toggleFunction=[self.ToggleFSIGFOrISIGIavailable, ['HKLIN_IS_I_SIGI', 'MERGED_OR_UNMERGED' ] ] )
    self.createLine( [ 'widget', '-browseDb', True, 'FREERFLAG' ] )
    self.createLine( [ 'widget', 'USE_TWIN', 'label', 'Twin refinement' ] )
    self.createLine( [ 'label', '<i>Warning: Intensities should be given for twin refinement. Using amplitudes is suboptimal.</i>' ],
                    toggleFunction=[self.ToggleTwinSuboptimal, ['HKLIN_IS_I_SIGI', 'F_SIGF_OR_I_SIGI', 'HKLIN']])
    self.closeSubFrame()

    self.openSubFrame(toggle = ['DATA_METHOD', 'open', [ 'spa' ] ] )
    self.createLine( [ 'label', 'Half map 1', 'widget', '-browseDb', True, 'MAPIN1' ] )
    self.createLine( [ 'label', 'Half map 2', 'widget', '-browseDb', True, 'MAPIN2' ] )
    self.createLine( [ 'label', 'Mask', 'widget', '-browseDb', True, 'MAPMASK' ] )
    self.createLine( [ 'label', 'Resolution:', 'stretch', 'widget', 'RES_MIN' ] )
    self.createLine( [ 'label', 'Mask radius:', 'stretch', 'widget', 'MASK_RADIUS' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Additional geometry dictionaries', 'stretch' ])
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', '-browseDb', True, 'DICT_LIST' ] )
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Options'] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'label', 'Number of refinement cycles:', 'stretch', 'widget', 'NCYCLES' ])#, toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']] )

    self.createLine( [ 'widget', 'HYDR_USE', 'label', 'Use riding hydrogens during refinement'], toggle = ['HYDR_USE', 'open', [ False ] ])
    self.createLine( [ 'widget', 'HYDR_USE', 'label', 'Use riding hydrogens during refinement', 'stretch', 'widget', 'HYDR_ALL'], toggle = ['HYDR_USE', 'open', [ True ] ] )
    add_waters = self.createLine( [ 'widget', 'ADD_WATERS', 'label', 'Add waters' ], toggle = ['DATA_METHOD', 'open', [ 'xtal' ] ])
    self.createLine( [ 'label', '&nbsp;and then perform further ', 'widget', 'NCYCLES_AFTER_ADD_WATERS', 'label', ' refinement cycles' ], toggle = [ 'ADD_WATERS','open', [ True ] ], appendLine=add_waters  )

    self.closeSubFrame()
    
    self.openFolder(folderFunction='controlParameters',title='Parameterisation', drawFolder=self.drawParameters)
    self.openFolder(folderFunction='controlParameters',title='Restraints', drawFolder=self.drawRestraints)
    self.openFolder(folderFunction='controlParameters',title='Advanced')
    self.drawAdvanced() # small change introduced to allow for automatically loading a keyword file in the 'advanced' tab
    self.setProsmartProteinMode()
    self.setProsmartNucleicAcidMode()
    self.setDnatcoMode()
    return

  def twinHelpPressed(self):
    msg = 'Twin refinement is not available for anomalous data.'

    text = 'To enable the twin refinement option,'
    text += ' select mean intensities or mean structure amplitudes'
    text += ' from the menu "Reflections".'
    text += '\n\n'
    text += 'Mean intensities are generated by\n'
    text += 'X-ray data reduction and analysis > Data reduction - AIMLESS\n'
    text += 'and are labeled as\n'
    text += '"... Mean intensities for twin refinement".'
    text += '\n\n'
    text += 'Mean structure amplitudes are generated by\n'
    text += 'Reflection Data Tools > Convert Intensities to Amplitudes\n'
    text += 'and are labeled as\n'
    text += '"... as Mean SFs".'

    msgBox = QtWidgets.QMessageBox(self)
    msgBox.setIcon(msgBox.Information)
    msgBox.setWindowModality(QtCore.Qt.WindowModal)
    msgBox.setText(msg)
    msgBox.setInformativeText(text)
    msgBox.setStandardButtons(msgBox.Ok)
    msgBox.setDefaultButton(msgBox.Ok)
    msgBox.exec_()

  def drawParameters( self ):
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
    self.createLine( [ 'subtitle', 'Atomic displacement parameters (ADPs)'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']] )
    self.createLine( [ 'widget', 'B_REFINEMENT_MODE', 'label', 'ADPs'] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Single particle analysis (SPA) settings'], toggle = ['DATA_METHOD', 'open', [ 'spa' ] ] )
    self.openSubFrame(frame=[True], toggle = ['DATA_METHOD', 'open', [ 'spa' ] ] )
    self.createLine( [ 'label', 'Pixel size (angstroem/pixel):', 'stretch', 'widget', 'PIXEL_SIZE' ] )
    self.createLine( [ 'label', 'Point group:', 'stretch', 'widget', 'POINTGROUP' ] )
    self.createLine( [ 'label', 'Helical twist (degrees):', 'stretch', 'widget', 'TWIST' ] )
    self.createLine( [ 'label', 'Helical rise (angstroem):', 'stretch', 'widget', 'RISE' ] )
    self.createLine( [ 'label', 'Set centre, <i>i.e.</i> the origin of symmetry. (Default is centre of the box.):', 'stretch', 'widget', 'CENTER_X', 'widget', 'CENTER_Y', 'widget', 'CENTER_Z' ] )
    self.createLine( [ 'label', 'Axis 1:', 'stretch', 'widget', 'AXIS1_X', 'widget', 'AXIS1_Y', 'widget', 'AXIS1_Z' ] )
    self.createLine( [ 'label', 'Axis 2:', 'stretch', 'widget', 'AXIS2_X', 'widget', 'AXIS2_Y', 'widget', 'AXIS2_Z' ] )
    self.createLine( [ 'label', '<i>Hint for axis: I: 5-fold, O: 4-fold, T: 3-fold, D<sub>n</sub>: 2-fold).</i>' ] )
    self.createLine( [ 'widget', 'IGNORE_SYMMETRY', 'label', 'Ignore symmetry information (MTRIX/_struct_ncs_oper) in the input structure model file' ])
    blurline = self.createLine( [ 'widget', 'BLURUSE', 'label', 'Map blurring (workaround for oversharpened maps)'] )
    self.createLine( [ 'label', indent + 'B-value:', 'stretch', 'widget', 'BLUR'], toggle = ['BLURUSE', 'open', [ True ] ] , appendLine=blurline)
    self.createLine( [ 'label', '<i>Note: This option does not affect output maps.</i>'], toggle = ['BLURUSE', 'open', [ True ] ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Scaling'], toggle = ['DATA_METHOD', 'open', [ 'xtal' ] ] )
    self.openSubFrame(frame=[True], toggle = ['DATA_METHOD', 'open', [ 'xtal' ] ] )
    self.createLine( [ 'widget', 'NO_SOLVENT', 'label', 'Do not consider bulk solvent contribution' ])
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Conformer groups and occupancy refinement'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']] )
    self.createLine( [ 'widget', 'OCCUPANCY_GROUPS', 'label', 'Specify partial occupancy groups (alternative conformers)' ] )
    self.createLine( [ 'widget', 'OCCUPANCY_SELECTION' ], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'OCCUPANCY_COMPLETE', 'label', 'Specify overlapping alternative conformer groups (constrain occupancies to sum to one)' ], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'OCCUPANCY_COMPLETE_TABLE' ], toggleFunction=[self.ToggleOccComplete, ['OCCUPANCY_GROUPS','OCCUPANCY_COMPLETE']] )
    self.createLine( [ 'widget', 'OCCUPANCY_INCOMPLETE', 'label', 'Specify overlapping alternative conformer groups (occupancies sum to less than one)' ], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'OCCUPANCY_INCOMPLETE_TABLE' ], toggleFunction=[self.ToggleOccIncomplete, ['OCCUPANCY_GROUPS','OCCUPANCY_INCOMPLETE']] )
    self.createLine( [ 'widget', 'OCCUPANCY_REFINEMENT', 'label', 'Perform refinement of atomic occupancies every', 'widget', 'OCCUPANCY_NCYCLE', 'label', 'cycle.'], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']] )
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()
    
    return

  @QtCore.Slot()
  def checkAllowUndefined(self):
    self.container.inputData.TLSIN.setQualifiers({'allowUndefined':str(self.container.controlParameters.TLSMODE) != 'FILE'})
    self.validate()
    return

  def drawRestraints( self ):
    self.createLine( [ 'subtitle', 'Weights'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    auto_weight = self.createLine( [ 'label', 'Weigh the experimental data using', 'widget', 'WEIGHT_OPT', 'label', 'weight'] )
    self.createLine( [ 'label', 'of', 'widget', 'WEIGHT' ], toggle = ['WEIGHT_OPT', 'open', [ 'MANUAL' ] ], appendLine=auto_weight )
    self.createLine( [ 'label', 'versus the restraints' ], appendLine=auto_weight )
    self.createLine( [ 'widget', 'WEIGHT_NO_ADJUST', 'label', 'Do not adjust weight during refinement'], toggle = ['WEIGHT_OPT', 'open', [ 'AUTO' ] ] )
    self.createLine( ['label', 'Bond RMSZ range for weight adjustment:', 'stretch',
                      'widget', 'WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN',
                      'widget', 'WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX'],
                      toggleFunction = [self.ToggleWeightAdjustRmszAvailable, ['WEIGHT_OPT', 'WEIGHT_NO_ADJUST', 'DATA_METHOD']])
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Non-Crystallographic Symmetry (NCS)'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOn,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    self.createLine( [ 'widget', 'USE_NCS', 'label', 'Use local non-crystallographic symmetry (NCS) restraints' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOff,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Covalent links'] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'FIND_LINKS', 'label', 'Detect and apply covalent linkages based on the current atomic coordinates' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Jelly-body'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleJellyOn,['UNRESTRAINED', 'FIX_XYZ']])
    use_jellybody = self.createLine( [ 'widget', 'USE_JELLY', 'label', 'Use jelly-body restraints' ] )
    self.createLine( [ 'label', '&nbsp;with sigma:', 'widget', 'JELLY_SIGMA', 'label', 'and max distance:', 'widget', 'JELLY_DIST' ], toggle = ['USE_JELLY', 'open', [ True ] ], appendLine=use_jellybody)
    self.createLine( [ 'widget', 'JELLY_ONLY', 'label', 'Jelly body refinement only' ], toggle = ['USE_JELLY', 'open', [ True ] ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleJellyOff,['UNRESTRAINED', 'FIX_XYZ']])
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'MetalCoord External Restraints for Metals'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOn,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    if shutil.which("metalCoord", mode=os.X_OK):
      if self.widget.subFrame is not None:
         self.currentLayout = self.widget.subFrame.layout()
      else:
         self.currentLayout = self.widget.currentFolderLayout
      self.ligands_checkboxes = ChoiceButtons()
      self.currentLayout.addWidget(self.ligands_checkboxes)
      self.updateMonomersWithMetalsWidget()
      self.createLine( [ 'widget', 'metalCoordPipeline.RUN_METALCOORD', 'label', 'Apply MetalCoord restraints for metal sites:' ] , toggle = ['metalCoordPipeline.RUN_METALCOORD', 'open', [ False ] ])
      self.createLine( [ 'widget', 'metalCoordPipeline.RUN_METALCOORD', 'label', 'Apply MetalCoord restraints for metal sites:' , 'stretch' , 'widget', 'metalCoordPipeline.GENERATE_OR_USE'] , toggle = ['metalCoordPipeline.RUN_METALCOORD', 'open', [ True ] ])
      self.createLine( [ 'widget', 'metalCoordPipeline.METALCOORD_RESTRAINTS'] , toggleFunction=[self.ToggleMetalCoordUse, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE']])
      self.ligands_checkboxes.clickedSignal.connect(self.updateMonomersWithMetalsSelection)
      self.createLine( [ 'widget', 'metalCoordPipeline.TOGGLE_ADVANCED', 'label', 'Show advanced options' ], toggleFunction=[self.ToggleMetalCoordGenerate, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE']] )
      self.createLine( [ 'label', 'Link records to metal sites in the atomic model:', 'stretch', 'widget', 'metalCoordPipeline.LINKS' ], toggle = ['metalCoordPipeline.TOGGLE_ADVANCED', 'open', [ True ] ] )
      self.createLine( [ 'label', 'Distance threshold: (range 0-1)<br/><i>A threshold d to select atoms is (r<sub>1</sub> + r<sub>2</sub>)*(1 + d) where r<sub>1</sub> and r<sub>2</sub> are covalent radii.</i>', 'stretch', 'widget', 'metalCoordWrapper.DISTANCE_THRESHOLD' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
      self.createLine( [ 'label', 'Maximum coordination number:', 'stretch', 'widget', 'metalCoordWrapper.MAXIMUM_COORDINATION_NUMBER' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
      self.createLine( [ 'label', 'Procrustes distance threshold: (range 0-1)', 'stretch', 'widget', 'metalCoordWrapper.PROCRUSTES_DISTANCE_THRESHOLD' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
      self.createLine( [ 'label', 'Minimum sample size for statistics:', 'stretch', 'widget', 'metalCoordWrapper.MINIMUM_SAMPLE_SIZE' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
      self.createLine( [ 'widget', 'metalCoordWrapper.USE_PDB', 'label', 'Use COD structures based on the input PDB/mmCIF coordinates' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
      self.createLine( [ 'widget', 'metalCoordWrapper.IDEAL_ANGLES', 'label', 'Provide only ideal bond angles' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
      self.createLine( [ 'widget', 'metalCoordWrapper.SIMPLE', 'label', 'Simple distance based filtering' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    else:
      self.createLine( [ 'label', '<i>MetalCoord is not installed.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOff,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'DNATCO External Restraints for Nucleic Acids'] )
    if self.isEditable():
       self.container.dnatco.TOGGLE_RESTRAINTS.dataChanged.connect(self.setDnatcoMode)
       self.container.controlParameters.REFINEMENT_MODE.dataChanged.connect(self.setDnatcoMode)
       self.container.controlParameters.UNRESTRAINED.dataChanged.connect(self.setDnatcoMode)
       self.container.controlParameters.FIX_XYZ.dataChanged.connect(self.setDnatcoMode)
       self.container.controlParameters.JELLY_ONLY.dataChanged.connect(self.setDnatcoMode)
    self.openSubFrame(frame=[True], toggle = ['dnatco.MODE', 'open', [ 'DISABLED' ] ] )
    self.createLine( [ 'label', '<i>Specify atomic model before setting up external restraints</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['dnatco.MODE', 'open', [ 'NONUCLEICACID' ] ] )
    self.createLine( [ 'label', '<i>Input atomic model contains no nucleotide chains</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['dnatco.MODE', 'open', [ 'RIGIDMODE' ] ] )
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['dnatco.MODE', 'open', [ 'NOTUSED' ] ] )
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['dnatco.MODE', 'open', [ 'UNSELECTED' ] ] )
    self.createLine( [ 'widget', 'dnatco.TOGGLE_RESTRAINTS', 'label', 'Generate and apply restraints for nucleic acids (experimental)' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['dnatco.MODE', 'open', [ 'SELECTED' ] ] )
    self.createLine( [ 'widget', 'dnatco.TOGGLE_RESTRAINTS', 'label', 'Generate and apply restraints for nucleic acids (experimental):' ] )
    self.autoGenerate(
         self.container.dnatco,
         selection={"includeParameters": ["MAX_RMSD", "RESTRAINTS_SIGMA"]},
      )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'ProSMART External Restraints for Protein Chains'] )
    if self.isEditable():
       self.container.prosmartProtein.TOGGLE.dataChanged.connect(self.setProsmartProteinMode)
       self.container.controlParameters.REFINEMENT_MODE.dataChanged.connect(self.setProsmartProteinMode)
       self.container.controlParameters.UNRESTRAINED.dataChanged.connect(self.setProsmartProteinMode)
       self.container.controlParameters.FIX_XYZ.dataChanged.connect(self.setProsmartProteinMode)
       self.container.controlParameters.JELLY_ONLY.dataChanged.connect(self.setProsmartProteinMode)
    self.openSubFrame(frame=[True], toggle = ['prosmartProtein.MODE', 'open', [ 'DISABLED' ] ] )
    self.createLine( [ 'label', '<i>Specify atomic model before setting up external restraints</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartProtein.MODE', 'open', [ 'NOPROTEIN' ] ] )
    self.createLine( [ 'label', '<i>Input atomic model contains no protein chains</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartProtein.MODE', 'open', [ 'RIGIDMODE' ] ] )
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartProtein.MODE', 'open', [ 'NOTUSED' ] ] )
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartProtein.MODE', 'open', [ 'UNSELECTED' ] ] )
    self.createLine( [ 'widget', 'prosmartProtein.TOGGLE', 'label', 'Generate and apply restraints for protein chains using homologous models' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartProtein.MODE', 'open', [ 'SELECTED' ] ] )
    self.createLine( [ 'widget', 'prosmartProtein.TOGGLE', 'label', 'Generate and apply restraints for protein chain(s):', 'widget', 'prosmartProtein.CHAINLIST_1', 'label', 'using homologous model(s):' ] )
    self.createLine( [ 'widget', '-browseDb', True, 'prosmartProtein.REFERENCE_MODELS' ] )
    self.createLine( [ 'label', 'Use', 'widget', 'prosmartProtein.ALL_BEST', 'label', 'chain(s) from the reference model(s).', 'stretch', 'label', 'Minimum sequence identity:', 'widget', 'prosmartProtein.SEQID', 'label', '%' ] )
    self.createLine( [ 'label', 'Generate restraints between', 'widget', 'prosmartProtein.SIDE_MAIN', 'label', 'atom-pairs.', 'stretch', 'label', 'Interatomic distance range:', 'widget', 'prosmartProtein.RMIN', 'label', 'to', 'widget', 'prosmartProtein.RMAX', 'label', 'angstroem' ] )
    self.createLine( [ 'label', 'Apply restraints up to distance', 'widget', 'prosmartProtein.DMAX', 'label', 'angstroem with robustness parameter (alpha)', 'widget', 'prosmartProtein.ALPHA' , 'stretch', 'label', 'Show advanced options', 'widget', 'prosmartProtein.ADVANCED' ])
    self.createLine( [ 'label', 'Minimum and maximum sigma', 'widget', 'prosmartProtein.SGMN', 'widget', 'prosmartProtein.SGMX'], toggle = ['prosmartProtein.ADVANCED', 'open', [ True ] ] )    
    self.createLine( [ 'widget', 'prosmartProtein.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high ADPs.' ] , toggleFunction=[self.hideProsmartProteinBfac, ['prosmartProtein.ADVANCED','prosmartProtein.TOGGLE_BFAC']])
    self.createLine( [ 'widget', 'prosmartProtein.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high ADPs.', 'stretch', 'label', 'Maximum B-factor: median plus ', 'widget', 'prosmartProtein.BFAC', 'label', 'x interquartile range' ] , toggleFunction=[self.showProsmartProteinBfac, ['prosmartProtein.ADVANCED','prosmartProtein.TOGGLE_BFAC']])
    self.createLine( [ 'widget', 'prosmartProtein.TOGGLE_ALT', 'label', 'Allow restraints involving atoms with alt codes.', 'stretch' , 'label', 'Ignore atoms with occupancies lower than', 'widget', 'prosmartProtein.OCCUPANCY'] , toggle = ['prosmartProtein.ADVANCED', 'open', [ True ] ] )
    self.createLine( [ 'label', 'Additional ProSMART keywords:', 'widget', 'prosmartProtein.KEYWORDS' ] , toggle = ['prosmartProtein.ADVANCED', 'open', [ True ] ] )
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'ProSMART External Restraints for Nucleic Acids'] )
    if self.isEditable():
       self.container.prosmartNucleicAcid.TOGGLE.dataChanged.connect( self.setProsmartNucleicAcidMode)
       self.container.controlParameters.REFINEMENT_MODE.dataChanged.connect(self.setProsmartNucleicAcidMode)
       self.container.controlParameters.UNRESTRAINED.dataChanged.connect(self.setProsmartNucleicAcidMode)
       self.container.controlParameters.FIX_XYZ.dataChanged.connect(self.setProsmartNucleicAcidMode)
       self.container.controlParameters.JELLY_ONLY.dataChanged.connect(self.setProsmartNucleicAcidMode)
    self.openSubFrame(frame=[True], toggle = ['prosmartNucleicAcid.MODE', 'open', [ 'DISABLED' ] ] )
    self.createLine( [ 'label', '<i>Specify atomic model before setting up external restraints</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartNucleicAcid.MODE', 'open', [ 'NONUCLEICACID' ] ] )
    self.createLine( [ 'label', '<i>Input atomic model contains no nucleotide chains</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartNucleicAcid.MODE', 'open', [ 'RIGIDMODE' ] ] )
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartNucleicAcid.MODE', 'open', [ 'NOTUSED' ] ] )
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartNucleicAcid.MODE', 'open', [ 'UNSELECTED' ] ] )
    self.createLine( [ 'widget', 'prosmartNucleicAcid.TOGGLE', 'label', 'Generate and apply restraints for nucleotide chain(s) using homologous models' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartNucleicAcid.MODE', 'open', [ 'SELECTED' ] ] )
    self.createLine( [ 'widget', 'prosmartNucleicAcid.TOGGLE', 'label', 'Generate and apply restraints for nucleotide chain(s):', 'widget', 'prosmartNucleicAcid.CHAINLIST_1', 'label', 'using homologous model(s):' ] )
    self.createLine( [ 'widget', '-browseDb', True, 'prosmartNucleicAcid.REFERENCE_MODELS' ] )
    self.createLine( [ 'label', 'Use', 'widget', 'prosmartNucleicAcid.ALL_BEST', 'label', 'chain(s) from the reference model(s).', 'stretch', 'label', 'Minimum sequence identity:', 'widget', 'prosmartNucleicAcid.SEQID', 'label', '%' ] )
    self.createLine( [ 'label', 'Generate restraints between', 'widget', 'prosmartNucleicAcid.SIDE_MAIN', 'label', 'atom-pairs.', 'stretch', 'label', 'Interatomic distance range:', 'widget', 'prosmartNucleicAcid.RMIN', 'label', 'to', 'widget', 'prosmartNucleicAcid.RMAX', 'label', 'angstroem' ] )
    self.createLine( [ 'label', 'Apply restraints up to distance', 'widget', 'prosmartNucleicAcid.DMAX', 'label', 'angstroem with weight', 'widget', 'prosmartNucleicAcid.WEIGHT', 'label', 'and robustness parameter (alpha)', 'widget', 'prosmartNucleicAcid.ALPHA', 'stretch', 'label', 'Show advanced options', 'widget', 'prosmartNucleicAcid.ADVANCED' ])
    self.createLine( [ 'widget', 'prosmartNucleicAcid.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high B-factors.' ], toggleFunction=[self.hideProsmartNucleicAcidBfac, ['prosmartNucleicAcid.ADVANCED','prosmartNucleicAcid.TOGGLE_BFAC']])
    self.createLine( [ 'widget', 'prosmartNucleicAcid.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high B-factors.', 'stretch', 'label', 'Maximum B-factor: median plus ', 'widget', 'prosmartNucleicAcid.BFAC', 'label', 'x interquartile range' ], toggleFunction=[self.showProsmartNucleicAcidBfac, ['prosmartNucleicAcid.ADVANCED','prosmartNucleicAcid.TOGGLE_BFAC']])
    self.createLine( [ 'widget', 'prosmartNucleicAcid.TOGGLE_ALT', 'label', 'Allow restraints involving atoms with alt codes.', 'stretch' , 'label', 'Ignore atoms with occupancies lower than', 'widget', 'prosmartNucleicAcid.OCCUPANCY'], toggle = ['prosmartNucleicAcid.ADVANCED', 'open', [ True ] ] )
    self.createLine( [ 'label', 'Additional ProSMART keywords:', 'widget','prosmartNucleicAcid.KEYWORDS' ], toggle = ['prosmartNucleicAcid.ADVANCED', 'open', [ True ] ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'ADP Restraints'] )
    self.openSubFrame()
    self.createLine( [ 'label', 'ADP restraint weight:', 'stretch', 'widget', 'ADPR_WEIGHT' ] )
    self.createLine( [ 'label', 'Maximum distance for ADP restraint:', 'stretch', 'widget', 'MAX_DIST_FOR_ADP_RESTRAINT' ] )
    self.createLine( [ 'label', 'ADP restraint power:', 'widget', 'stretch', 'ADP_RESTRAINT_POWER' ] )
    self.createLine( [ 'widget', 'ADP_RESTRAINT_NO_LONG_RANGE', 'label', 'No long range for ADP restraint' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Van der Waals Repulsion Restraints'] )
    self.openSubFrame()
    self.createLine( [ 'label', 'Van der Waals repulsion restraint weight:', 'stretch', 'widget', 'VDWR_WEIGHT' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Infrequently Used Options'] )
    self.openSubFrame()
    self.createLine( [ 'widget', 'UNRESTRAINED', 'label', 'No positional restraints' ] )
    self.createLine( [ 'widget', 'FIX_XYZ', 'label', 'Fix coordinates' ] )
    self.closeSubFrame()
    return

  def showProsmartProteinBfac(self):
     return bool(self.container.prosmartProtein.TOGGLE_BFAC)

  def showProsmartNucleicAcidBfac(self):
     return bool(self.container.prosmartNucleicAcid.TOGGLE_BFAC)

  def hideProsmartProteinBfac(self):
     return self.container.prosmartProtein.ADVANCED and not self.container.prosmartProtein.TOGGLE_BFAC

  def hideProsmartNucleicAcidBfac(self):
     return self.container.prosmartNucleicAcid.ADVANCED and not self.container.prosmartNucleicAcid.TOGGLE_BFAC

  def drawAdvanced( self ):
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

    self.openSubFrame(frame=[True], toggle = ['DATA_METHOD', 'open', [ 'xtal' ] ] )
    custom_res = self.createLine( [ 'widget', 'RES_CUSTOM', 'label', 'Use custom resolution limits' ] )
    self.createLine( [ 'label', indent+indent+'highest (d<sub>min</sub>):', 'widget', 'RES_MIN', 'label', ' lowest (d<sub>max</sub>):', 'widget', 'RES_MAX' ], toggle = ['RES_CUSTOM', 'open', [ True ] ], appendLine=custom_res )
    self.createLine( [ 'label', 'FreeR flag number for test set:' , 'widget', 'FREERFLAG_NUMBER'])
    self.createLine( [ 'label', 'Diffraction experiment type:', 'widget', 'SCATTERING_FACTORS' ] )
    if self.isEditable():
        self.container.controlParameters.SCATTERING_FACTORS.dataChanged.connect(self.ExperimentChanged)
    self.createLine( [ 'widget', 'USE_WORK_IN_EST', 'label', 'Use work reflections in maximum likelihood parameter estimates' ] )
    self.closeSubFrame()

    self.createLine( [ 'widget', 'CROSS_VALIDATION', 'label', 'Cross validation with half maps' ], toggle = ['DATA_METHOD', 'open', [ 'spa' ] ] )

    self.createLine( [ 'widget', 'H_OUT', 'label', 'Write hydrogen atoms in the output model' ] )
    self.createLine( [ 'widget', 'H_REFINE', 'label', 'Refine hydrogen positions' ], toggle = ['HYDR_USE', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'KEEP_CHARGES', 'label', 'Keep charges, i.e. use scattering factor for charged atoms where relevant' ] )

    self.createLine( [ 'subtitle', 'Structure model modification before refinement' ] )
    self.openSubFrame(frame=[True])
    reset_bfac = self.createLine( [ 'widget', 'BFACSETUSE', 'label', 'Reset all ADPs at start' ])
    self.createLine( [ 'label', '&nbsp;to fixed value:', 'widget', 'BFACSET' ], toggle = ['BFACSETUSE', 'open', [ True ] ], appendLine=reset_bfac )
    randomize = self.createLine( [ 'widget', 'RANDOMIZEUSE', 'label', 'Shake coordinates at start' ])
    self.createLine( [ 'label', '&nbsp;with specified RMSD:', 'widget', 'RANDOMIZE' ], toggle = ['RANDOMIZEUSE', 'open', [ True ] ], appendLine=randomize )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Additional keywords'] )
    self.createLine( [ 'widget', '-browseDb', True, 'SERVALCAT_KEYWORD_FILE' ] )
    self.createLine( [ 'label', 'Extra servalcat command line options:', 'widget', 'EXTRA_SERVALCAT_OPTIONS' ] )
    self.getWidget('EXTRA_SERVALCAT_OPTIONS').setFixedWidth(400)

    self.createLine( [ 'subtitle', 'Validation and Analysis' ] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'VALIDATE_IRIS', 'label', 'Generate Iris validation report' ] )
    self.createLine( [ 'widget', 'VALIDATE_RAMACHANDRAN', 'label', 'Generate Ramachandran plots' ] )
    self.createLine( [ 'widget', 'VALIDATE_MOLPROBITY', 'label', 'Run MolProbity to analyse geometry' ] )
    self.createLine( [ 'widget', 'dnatco.TOGGLE_VALIDATION', 'label', 'Run DNATCO to validate nucleic acids' ] )

    self.createLine( [ 'widget', 'RUN_ADP_ANALYSIS', 'label', 'Run ADP analysis' ] )
    self.createLine( [ 'label', 'Atoms with a B-value lower than <i>the first quartile - factor * interquartile_range</i><br />or higher than <i>the third quartile + factor * interquartile_range</i> to be reported. Factor:',
                       'stretch', 'widget', 'ADP_IQR_FACTOR' ], toggle = ['RUN_ADP_ANALYSIS', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'RUN_COORDADPDEV_ANALYSIS', 'label', 'Run analysis of changes in coordinates and ADPs' ] )
    self.createLine( [ 'label', 'Minimum shift of atom coordinates to be reported:', 'stretch', 'widget', 'monitor.MIN_COORDDEV' ], toggle = ['RUN_COORDADPDEV_ANALYSIS', 'open', [ True ] ] )
    self.createLine( [ 'label', 'Minimum shift of B-values to be reported:', 'stretch', 'widget', 'monitor.MIN_ADPDEV' ], toggle = ['RUN_COORDADPDEV_ANALYSIS', 'open', [ True ] ] )
    self.closeSubFrame()
    return

#-  --------------------          --------------------          --------------------

  def anomalousAvailable(self):
    if not self.container.inputData.HKLIN.isSet(): return False
    if not self.isEditable():
       #only display after job is run if USEANOMALOUS has already been set to True
       return bool(self.container.controlParameters.USEANOMALOUS)
    #Peak to see if we can make F+/F-
    self.container.inputData.HKLIN.setContentFlag(reset=True)
    canConvertString, toType = self.container.inputData.HKLIN.conversion(2)
    print('Can convert ?',self.container.inputData.HKLIN.contentFlag, canConvertString, toType)
    self.validate()
    if canConvertString == 'no':
       self.container.controlParameters.USEANOMALOUS = False
       return False
    else:
       self.container.controlParameters.USEANOMALOUS = True
       return True

  @QtCore.Slot()
  def hklinChanged(self):
    self.getObsType()

  @QtCore.Slot()
  def ExperimentChanged(self):
    if str(self.container.controlParameters.SCATTERING_FACTORS) == "NEUTRON":
      self.container.controlParameters.H_REFINE.set(True)
    else:
      self.container.controlParameters.H_REFINE.set(False)

  @QtCore.Slot()
  def modelChanged(self):
     self.container.prosmartProtein.TOGGLE.set(False)
     self.container.prosmartProtein.CHAINLIST_1.unSet()
     self.container.prosmartNucleicAcid.TOGGLE.set(False)
     self.container.prosmartNucleicAcid.CHAINLIST_1.unSet()
     self.container.inputData.AMINOACID_CHAINS.unSet()
     self.container.inputData.NUCLEOTIDE_CHAINS.unSet()
     if self.container.inputData.XYZIN.isSet():
        chain_list = self.getChainList()
        if chain_list:
           self.container.inputData.AMINOACID_CHAINS.set(chain_list["aminoacid"])
           self.container.inputData.NUCLEOTIDE_CHAINS.set(chain_list["nucleotide"])
           chain_list_protein = ' '.join(chain_list["aminoacid"])
           chain_list_nucleicacid = ' '.join(chain_list["nucleotide"])
           self.container.prosmartProtein.CHAINLIST_1.set(chain_list_protein)
           self.container.prosmartNucleicAcid.CHAINLIST_1.set(chain_list_nucleicacid)
           self.updateViewFromModel()
     self.setProsmartProteinMode()
     self.setProsmartNucleicAcidMode()
     self.setDnatcoMode()
     self.getMonomersWithMetals()

  @QtCore.Slot()
  def setProsmartProteinMode(self):
     self.container.prosmartProtein.CHAINLIST_1.setQualifiers({'allowUndefined':True})
     self.container.prosmartProtein.REFERENCE_MODELS.setQualifiers({'listMinLength':0})
     if self.ToggleRigidModeOn():
        self.container.prosmartProtein.MODE.set('RIGIDMODE')
     elif self.ToggleRestraintsOff():
        self.container.prosmartProtein.MODE.set('NOTUSED')
     elif self.container.inputData.XYZIN.isSet():
        if self.container.prosmartProtein.CHAINLIST_1.isSet():
           if len(self.container.prosmartProtein.CHAINLIST_1.__str__())>0:
              if self.container.prosmartProtein.TOGGLE:
                 self.container.prosmartProtein.MODE.set('SELECTED')
                 self.container.prosmartProtein.CHAINLIST_1.setQualifiers({'allowUndefined':False})
                 self.container.prosmartProtein.REFERENCE_MODELS.setQualifiers({'listMinLength':1})
              else:
                 self.container.prosmartProtein.MODE.set('UNSELECTED')
              self.validate()
              return
        self.container.prosmartProtein.MODE.set('NOPROTEIN')
     else:
        self.container.prosmartProtein.MODE.set('DISABLED')
     self.validate()

  @QtCore.Slot()
  def setProsmartNucleicAcidMode(self):
     self.container.prosmartNucleicAcid.CHAINLIST_1.setQualifiers({'allowUndefined':True})
     self.container.prosmartNucleicAcid.REFERENCE_MODELS.setQualifiers({'listMinLength':0})
     if self.ToggleRigidModeOn():
        self.container.prosmartNucleicAcid.MODE.set('RIGIDMODE')
     elif self.ToggleRestraintsOff():
        self.container.prosmartNucleicAcid.MODE.set('NOTUSED')
     elif self.container.inputData.XYZIN.isSet():
        if self.container.prosmartNucleicAcid.CHAINLIST_1.isSet():
           if len(self.container.prosmartNucleicAcid.CHAINLIST_1.__str__())>0:
              if self.container.prosmartNucleicAcid.TOGGLE:
                 self.container.prosmartNucleicAcid.MODE.set('SELECTED')
                 self.container.prosmartNucleicAcid.CHAINLIST_1.setQualifiers({'allowUndefined':False})
                 self.container.prosmartNucleicAcid.REFERENCE_MODELS.setQualifiers({'listMinLength':1})
              else:
                 self.container.prosmartNucleicAcid.MODE.set('UNSELECTED')
              self.validate()
              return
        self.container.prosmartNucleicAcid.MODE.set('NONUCLEICACID')
     else:
        self.container.prosmartNucleicAcid.MODE.set('DISABLED')
     self.validate()

  @QtCore.Slot()
  def setDnatcoMode(self):
     if self.ToggleRigidModeOn():
        self.container.dnatco.MODE.set('RIGIDMODE')
     elif self.container.inputData.XYZIN.isSet():
        if self.container.inputData.NUCLEOTIDE_CHAINS.isSet():
           if self.container.inputData.NUCLEOTIDE_CHAINS:
              if self.container.dnatco.TOGGLE_RESTRAINTS:
                 self.container.dnatco.MODE.set('SELECTED')
              else:
                 self.container.dnatco.MODE.set('UNSELECTED')
              self.container.dnatco.TOGGLE_VALIDATION.set(True)
              self.validate()
              return
        self.container.dnatco.MODE.set('NONUCLEICACID')
        self.container.dnatco.TOGGLE_VALIDATION.set(False)
     else:
        self.container.dnatco.MODE.set('DISABLED')
     self.validate()

  def getMonomersWithMetals(self):
        if self.container.inputData.XYZIN.isSet():
            if os.path.isfile(str(self.container.inputData.XYZIN.fullPath)):
                self.monomersWithMetals = set()
                try:
                    st = gemmi.read_structure(str(self.container.inputData.XYZIN.fullPath))
                    for model in st:
                        lookup = {x.atom: x for x in model.all()}
                        for cra in lookup.values():
                            if cra.atom.element.is_metal:
                                if cra.residue.name not in self.monomersWithMetals:
                                    self.monomersWithMetals.add(cra.residue.name)
                    self.monomersWithMetals = list(self.monomersWithMetals)
                except Exception as e:
                    print("Getting codes for monomers containing metals was not successful: ", e)
                    self.monomersWithMetals = []
                self.container.metalCoordPipeline.LIGAND_CODES_AVAILABLE.set(self.monomersWithMetals)
                self.updateMonomersWithMetalsWidget()

  def updateMonomersWithMetalsWidget(self):
        if hasattr(self, "ligands_checkboxes"):
            if self.monomersWithMetals:
                self.ligands_checkboxes.setChoices(
                    "Codes of monomers including metal sites:",
                    self.monomersWithMetals, tags=[], notes=[], exclusiveChoice=False)
            else:
                self.ligands_checkboxes.setChoices(
                    "<i>No monomers including metal sites were found in the input atomic model.</i>",
                    self.monomersWithMetals, tags=[], notes=[], exclusiveChoice=False)

  def updateMonomersWithMetalsSelection(self):
        if hasattr(self, "ligands_checkboxes"):   \
            self.container.metalCoordPipeline.LIGAND_CODES_SELECTED = self.ligands_checkboxes.selectedList

  def getChainList(self):
     chain_list = {}
     if self.container.inputData.XYZIN.isSet():
         model_path = self.container.inputData.XYZIN.fullPath.__str__()
         import mmdb2 as mmdb
         molHnd = mmdb.Manager()
         molHnd.SetFlag(mmdb.MMDBF_IgnoreRemarks)
         molHnd.SetFlag(mmdb.MMDBF_IgnoreBlankLines)
         molHnd.SetFlag(mmdb.MMDBF_IgnoreHash)
         molHnd.SetFlag(mmdb.MMDBF_IgnoreNonCoorPDBErrors)
         molHnd.ReadCoorFile(model_path)
         if molHnd.GetNumberOfModels() > 0:
            chain_list = {"aminoacid":[],"nucleotide":[],"other":[]}
            model = molHnd.GetModel(1)
            for chainIdx in range(0,model.GetNumberOfChains()):
               chain = model.GetChain(chainIdx)
               if chain.isAminoacidChain():
                  chain_list["aminoacid"].append(chain.GetChainID())
               elif chain.isNucleotideChain():
                  chain_list["nucleotide"].append(chain.GetChainID())
               else:
                  chain_list["other"].append(chain.GetChainID())
     return chain_list

  def getObsType(self):
    if self.container.inputData.HKLIN.isSet():
      if self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR:
        self.container.controlParameters.HKLIN_IS_I_SIGI = True
      elif self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN:
        self.container.controlParameters.HKLIN_IS_I_SIGI = True
      else:
        self.container.controlParameters.HKLIN_IS_I_SIGI = False
        self.container.controlParameters.F_SIGF_OR_I_SIGI = 'F_SIGF'

  def isValid(self):
      invalidElements = super(Cservalcat_pipe, self).isValid()
      #Check whether invocation is from runTask
      import traceback
      functionNames = [a[2] for a in traceback.extract_stack()]

      if str(self.container.controlParameters.DATA_METHOD) == 'xtal':
         if self.container.inputData.HKLIN.isSet() and  self.container.inputData.XYZIN.isSet() and self.container.inputData.XYZIN.fileContent.mmdbManager and hasattr(self.container.inputData.XYZIN.fileContent.mmdbManager,"GetCell"):
            cellMismatch = False
            sgMismatch = False
            if abs(float(self.container.inputData.HKLIN.fileContent.cell.a) - self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()[1]) > 1e-2: cellMismatch = True
            if abs(float(self.container.inputData.HKLIN.fileContent.cell.b) - self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()[2]) > 1e-2: cellMismatch = True
            if abs(float(self.container.inputData.HKLIN.fileContent.cell.c) - self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()[3]) > 1e-2: cellMismatch = True
            if abs(float(self.container.inputData.HKLIN.fileContent.cell.alpha) - self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()[4]) > 1e-2: cellMismatch = True
            if abs(float(self.container.inputData.HKLIN.fileContent.cell.beta) - self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()[5]) > 1e-2: cellMismatch = True
            if abs(float(self.container.inputData.HKLIN.fileContent.cell.gamma) - self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()[6]) > 1e-2: cellMismatch = True
            if str(self.container.inputData.XYZIN.fileContent.mmdbManager.GetSpaceGroup()) != str(self.container.inputData.HKLIN.fileContent.spaceGroup): sgMismatch = True
            if cellMismatch or sgMismatch:
               if functionNames[-2] == 'runTask':
                     from PySide2.QtWidgets import QMessageBox
                     msg = QMessageBox()
                     msg.setIcon(QMessageBox.Question)
                     msg.setText("Warning")
                     infoText = "You are trying to execute refinement with possible space group/cell mismatch between reflection and model.<br/>"
                     infoText += "<br/>Reflections:<br/>"
                     refCell = self.container.inputData.HKLIN.fileContent.cell
                     infoText += "%.3f %.3f %.3f<br/>%.3f %.3f %.3f<br/>%s<br/>" % (float(refCell.a),  float(refCell.b), float(refCell.c) ,float(refCell.alpha), float(refCell.beta), float(refCell.gamma),  str(self.container.inputData.HKLIN.fileContent.spaceGroup))
                     infoText += "<br/>Model:<br/>"
                     modCell = self.container.inputData.XYZIN.fileContent.mmdbManager.GetCell()
                     infoText += "%.3f %.3f %.3f<br/>%.3f %.3f %.3f<br/>%s<br/>" % (float(modCell[1]), float(modCell[2]), float(modCell[3]),float(modCell[4]), float(modCell[5]), float(modCell[6]), str(self.container.inputData.XYZIN.fileContent.mmdbManager.GetSpaceGroup()))
                     msg.setInformativeText(infoText)
                     msg.setStandardButtons(QMessageBox.Ignore | QMessageBox.Cancel)
                     retval = msg.exec_()
                     if retval == QMessageBox.Cancel:
                        invalidElements.append(self.container.inputData.HKLIN)
                        invalidElements.append(self.container.inputData.XYZIN)

         if functionNames[-2] == 'runTask':
            if not self.container.inputData.FREERFLAG.isSet():
               from PySide2.QtWidgets import QMessageBox
               msg = QMessageBox()
               msg.setIcon(QMessageBox.Question)

               msg.setText("Warning")
               msg.setInformativeText("You are trying to execute refinement without selecting a Free-R data item")
               msg.setWindowTitle("Free-R warning")
               msg.setDetailedText("While refinement without FreeR is reasonable under some circumstances, it is generally discouraged because it can provide a misleading impression of progress in refinement where over-fitting of the data can occur")
               msg.setStandardButtons(QMessageBox.Ignore | QMessageBox.Cancel)
               retval = msg.exec_()
               if retval == QMessageBox.Cancel:
                  invalidElements.append(self.container.inputData.FREERFLAG)

         if self.container.controlParameters.REFINEMENT_MODE.isSet():
            if str(self.container.controlParameters.REFINEMENT_MODE) == 'RIGID':
               for sel0 in self.container.controlParameters.RIGID_BODY_SELECTION:
                  sel = sel0.get()
                  if not sel['groupId'] or not sel['chainId'] or not sel['firstRes'] or not sel['lastRes']:
                     invalidElements.append(self.container.controlParameters.RIGID_BODY_SELECTION)
            else:
               if self.container.controlParameters.OCCUPANCY_GROUPS:
                  occup_groupids = []
                  for sel0 in self.container.controlParameters.OCCUPANCY_SELECTION:
                     sel = sel0.get()
                     if not sel['groupId'] or not sel['chainIds'] or not sel['firstRes'] or not sel['lastRes']:
                        invalidElements.append(self.container.controlParameters.OCCUPANCY_SELECTION)
                     if sel['groupId']:
                        occup_groupids.append(str(sel['groupId']))
                  if self.container.controlParameters.OCCUPANCY_COMPLETE:
                     for sel0 in self.container.controlParameters.OCCUPANCY_COMPLETE_TABLE:
                        sel = sel0.get()
                        for id in sel['groupIds'].split(' '):
                           if not id in occup_groupids:
                              invalidElements.append(self.container.controlParameters.OCCUPANCY_COMPLETE_TABLE)
                  if self.container.controlParameters.OCCUPANCY_INCOMPLETE:
                     for sel0 in self.container.controlParameters.OCCUPANCY_INCOMPLETE_TABLE:
                        sel = sel0.get()
                        for id in sel['groupIds'].split(' '):
                           if not id in occup_groupids:
                              invalidElements.append(self.container.controlParameters.OCCUPANCY_INCOMPLETE_TABLE)

      elif str(self.container.controlParameters.DATA_METHOD) == 'spa':
         if not self.container.inputData.MAPIN1.isSet() or not self.container.inputData.MAPIN2.isSet():
            if functionNames[-2] == 'runTask':
               from PySide2.QtWidgets import QMessageBox
               msg = QMessageBox()
               msg.setIcon(QMessageBox.Critical)
               msg.setText("Error")
               msg.setInformativeText("Two half maps are required but were not provided.")
               msg.setWindowTitle("Half maps missing")
               msg.setStandardButtons(QMessageBox.Cancel)
               retval = msg.exec_()
               if not self.container.inputData.MAPIN1.isSet():
                  invalidElements.append(self.container.inputData.MAPIN1)
               elif not self.container.inputData.MAPIN2.isSet():
                  invalidElements.append(self.container.inputData.MAPIN2)
         if not self.container.inputData.MAPMASK.isSet():
            if functionNames[-2] == 'runTask':
               from PySide2.QtWidgets import QMessageBox
               msg = QMessageBox()
               msg.setIcon(QMessageBox.Critical)
               msg.setText("Error")
               msg.setInformativeText("Map mask is required but was not provided.")
               msg.setWindowTitle("Map mask missing")
               msg.setStandardButtons(QMessageBox.Cancel)
               retval = msg.exec_()
               invalidElements.append(self.container.inputData.MAPMASK)
         if not self.container.controlParameters.RES_MIN.isSet():
            if functionNames[-2] == 'runTask':
               from PySide2.QtWidgets import QMessageBox
               msg = QMessageBox()
               msg.setIcon(QMessageBox.Critical)
               msg.setText("Error")
               msg.setInformativeText("Resolution is required but was not provided.")
               msg.setWindowTitle("Resolution missing")
               msg.setStandardButtons(QMessageBox.Cancel)
               retval = msg.exec_()
               invalidElements.append(self.container.controlParameters.RES_MIN)

      return invalidElements

