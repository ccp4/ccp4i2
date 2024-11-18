from __future__ import print_function

#from lxml import etree

"""
    servalcat_xtal_pipe_gui.py: CCP4 GUI Project
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

"""
     Andrey Lebedev September 2011 - refmac_martin gui
     Liz Potterton Aug 2012 - convert for MTZ ADO's demo
     Liz Potterton Oct 2012 - Moved mini-MTZ version to refmac_martin
"""

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets
from qtgui import CCP4XtalWidgets
from core.CCP4PluginScript import CPluginScript
from core import CCP4XtalData
from pipelines.import_merged.script.dybuttons import *  # for ChoiceButtons()
import os
import gemmi


def whatNext(jobId=None,childTaskName=None,childJobNumber=None,projectName=None):
    import os
    from core import CCP4Modules, CCP4Utils, CCP4File, CCP4Container, CCP4Data, CCP4PluginScript
    jobStatus = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId,'status')
    if jobStatus == 'Unsatisfactory':
        returnList = ['LidiaAcedrg', 'servalcat_xtal_pipe']
    else:
        returnList = ['servalcat_xtal_pipe', 'coot_rebuild', 'moorhen_rebuild', 'modelcraft']
    return returnList

class Cservalcat_xtal_pipe(CCP4TaskWidget.CTaskWidget):

  TASKTITLE='Refinement - Servalcat (experimental)'
  SHORTTASKTITLE='servalcat_xtal_pipe'
  # DESCRIPTION='Refine (Refmacat/Refmac5) with optional restraints (Prosmart, Platonyzer)'
  TASKNAME = 'servalcat_xtal_pipe'
  TASKLABEL = 'servalcat_xtal_pipe'
  TASKVERSION = 0.0
  TASKMODULE = 'refinement'
  MGDISPLAYFILES = ['XYZOUT','FPHIOUT','DIFFPHIOUT']
  AUTOPOPULATEINPUT = True
  RANK=1

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def ToggleWeightAuto(self):
    return str(self.container.controlParameters.WEIGHT_OPT) == 'MANUAL'

  #def ToggleTLS(self):
  #return self.container.inputData.TLSIN.isSet()
  #def ToggleTLSNot(self):
  #return not self.container.inputData.TLSIN.isSet()
  #def ToggleTLSUsed(self):
  #return self.container.inputData.TLSIN.isSet() or self.container.controlParameters.AUTOTLS
  #def ToggleTLSNotUsed(self):
  #return not self.container.inputData.TLSIN.isSet() and not self.container.controlParameters.AUTOTLS

  def ToggleTLSModeOn(self):
    if str(self.container.controlParameters.TLSMODE) != 'NONE':
        from core import CCP4Modules
        CurrentStatus = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(self._jobId,'status')
        if CurrentStatus == "Pending":
            self.container.controlParameters.BFACSETUSE = True
        return True
    return False

  def ToggleTLSModeOff(self):
    if str(self.container.controlParameters.TLSMODE) == 'NONE':
        from core import CCP4Modules
        CurrentStatus = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(self._jobId,'status')
        if CurrentStatus == "Pending":
            self.container.controlParameters.BFACSETUSE = False
        return True
    return False

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

  def ToggleNeutronModeOn(self):
    return str(self.container.controlParameters.SCATTERING_FACTORS) == 'NEUTRON'

  def ToggleNeutronModeOff(self):
    return not str(self.container.controlParameters.SCATTERING_FACTORS) == 'NEUTRON'
    
  def ToggleNeutronModeHD(self):
    if not self.container.controlParameters.HYDR_USE:
      return False
    return self.container.controlParameters.HD_FRACTION

  def ToggleNeutronModeHD_YES(self):
    if not self.container.controlParameters.HYDR_USE:
      return False
    if str(self.container.controlParameters.HYDR_ALL) == "ALL":
      return False
    return self.container.controlParameters.HD_FRACTION

  def ToggleNeutronModeHD_ALL(self):
    if not self.container.controlParameters.HYDR_USE:
      return False
    if str(self.container.controlParameters.HYDR_ALL) == "YES":
      return False
    return self.container.controlParameters.HD_FRACTION

  def ToggleH_REFINE(self):
    if not self.container.controlParameters.HYDR_USE:
      return False
    return self.container.controlParameters.H_REFINE

  def ToggleTLSModeFile(self):
    return str(self.container.controlParameters.TLSMODE) == 'FILE'

  #def ToggleTwinAvailable(self):
  #  if not self.container.controlParameters.USEANOMALOUS:
  #      return True
  #  self.container.controlParameters.USE_TWIN = False
  #  return False

  #def ToggleTwinNotAvailable(self):
  # if self.container.controlParameters.USEANOMALOUS:
  #    return True
  # return False

  def ToggleTwinSuboptimal(self):
    return (not self.container.controlParameters.HKLIN_IS_I_SIGI or self.container.controlParameters.F_SIGF_OR_I_SIGI == "F_SIGF")

  def CheckScaleType(self):
    if str(self.container.controlParameters.SCALE_TYPE) == 'BABINET':
        self.container.controlParameters.SOLVENT_MASK_TYPE.set('NO')
    else:
        self.container.controlParameters.SOLVENT_MASK_TYPE.set('EXPLICIT')
    return True

  def ToggleUseSolventMask(self):
    return str(self.container.controlParameters.SOLVENT_MASK_TYPE) == 'EXPLICIT'

  def ToggleSolventAdvanced(self):
    return str(self.container.controlParameters.SOLVENT_MASK_TYPE) == 'EXPLICIT' and self.container.controlParameters.SOLVENT_ADVANCED

  def ToggleElectronDiffraction(self):
    return str(self.container.controlParameters.SCATTERING_FACTORS) == 'ELECTRON'

  def ToggleOccComplete(self):
    if self.container.controlParameters.OCCUPANCY_GROUPS and self.container.controlParameters.OCCUPANCY_COMPLETE:
        return True
    return False

  def ToggleOccIncomplete(self):
    if self.container.controlParameters.OCCUPANCY_GROUPS and self.container.controlParameters.OCCUPANCY_INCOMPLETE:
        return True
    return False

  def TogglePlatonyzerNAMG(self):
    if self.container.platonyzer.TOGGLE and str(self.container.platonyzer.MODE) == 'NA_MG':
        return True
    return False

  def drawContents(self):
    self.setProgramHelpFile('servalcat_xtal')
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
    if self.container.metalCoordPipeline.LIGAND_CODES_AVAILABLE:
        self.monomersWithMetals = self.container.metalCoordPipeline.LIGAND_CODES_AVAILABLE
        self.container.metalCoordPipeline.LIGAND_CODES_SELECTED = self.container.metalCoordPipeline.LIGAND_CODES_AVAILABLE
    else:
        self.monomersWithMetals = []
    #-  --------------------          --------------------          --------------------
    folder = self.openFolder(folderFunction='inputData',title='Input Data')
    self.hklinChanged()

    #self.createLine( [ 'advice',' '] )
    self.createLine( [ 'subtitle', 'Main inputs' ])
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', '-browseDb', True, 'XYZIN' ])
    self.closeSubFrame()
    if self.isEditable():
       self.container.inputData.XYZIN.dataChanged.connect( self.modelChanged)
       #try:
       #self.getWidget('XYZIN').showAtomSelection()
       #except:
       #pass
    self.createLine( [ 'label', 'Experimental data type:', 'widget', 'DATA_METHOD' ])
    self.openSubFrame(toggle = ['DATA_METHOD', 'open', [ 'xtal' ] ] )
    self.createLine( [ 'widget', '-browseDb', True, 'HKLIN' ] )
    self.container.inputData.HKLIN.dataChanged.connect( self.hklinChanged )
    self.createLine( [ 'label', 'Refinement against <b>amplitudes</b>.'], toggle = ['HKLIN_IS_I_SIGI', 'open', [ False ] ] )
    self.createLine( [ 'label', 'Refinement against', 'widget', 'F_SIGF_OR_I_SIGI'], toggle = ['HKLIN_IS_I_SIGI', 'open', [ True ] ] )
    #self.createLine( [ 'widget', '-browseDb', True, 'UNMERGED' ])
    #self.createLine( [ 'label','Use anomalous data for ', 'widget', 'USEANOMALOUSFOR', 'stretch', 'label', 'Wavelength', 'widget', 'WAVELENGTH'],toggleFunction=[self.anomalousAvailable,['HKLIN']])
    #if self.isEditable():
    #    if not self.container.controlParameters.WAVELENGTH.isSet(): self.getWavelength()
    self.createLine( [ 'widget', '-browseDb', True, 'FREERFLAG' ] )
    self.closeSubFrame()
    self.openSubFrame(toggle = ['DATA_METHOD', 'open', [ 'spa' ] ] )
    self.createLine( [ 'label', 'Half map 1', 'widget', '-browseDb', True, 'MAPIN1' ] )
    self.createLine( [ 'label', 'Half map 2', 'widget', '-browseDb', True, 'MAPIN2' ] )
    self.createLine( [ 'label', 'Mask', 'widget', '-browseDb', True, 'MAPMASK' ] )
    self.createLine( [ 'label', 'Resolution', 'widget', 'RES_MIN' ] )
  
    self.closeSubFrame()
    #self.createLine( [ 'advice',' '] )
    #self.createLine( [ 'subtitle', 'Experimental phase information', 'stretch' ])
    #self.openSubFrame(frame=[True])
    #self.createLine( [ 'widget', '-browseDb', True, 'ABCD' ] )
    #self.closeSubFrame()
    self.createLine( [ 'subtitle', 'Additional geometry dictionaries', 'stretch' ])
    self.openSubFrame(frame=[True])
    #self.createLine( [ 'widget', '-browseDb', True, 'DICT' ] )
    self.createLine( [ 'widget', '-browseDb', True, 'DICT_LIST' ] )
    #if self.isEditable():
       #self.container.inputData.DICT_LIST.dataChanged.connect( self.MergeDictionaries)
    self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Options'] )
    self.openSubFrame(frame=[True])
    #self.createLine( [ 'label', 'Refinement mode (restrained or rigid body):', 'stretch', 'widget', 'REFINEMENT_MODE' ] )
    self.createLine( [ 'label', 'Number of refinement cycles:', 'stretch', 'widget', 'NCYCLES' ])#, toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']] )
    #self.createLine( [ 'label', 'Number of rigid body refinement cycles:', 'stretch', 'widget', 'NCYCRIGID' ], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']] )

    #use_hydr = self.createLine( [ 'widget', 'HYDR_USE', 'label', 'Use riding hydrogens during refinement'] )
    #self.createLine( [ 'stretch', 'widget', 'HYDR_ALL'], toggle = ['HYDR_USE', 'open', [ True ] ] , appendLine=use_hydr)
    self.createLine( [ 'widget', 'HYDR_USE', 'label', 'Use riding hydrogens during refinement'], toggle = ['HYDR_USE', 'open', [ False ] ])
    self.createLine( [ 'widget', 'HYDR_USE', 'label', 'Use riding hydrogens during refinement', 'stretch', 'widget', 'HYDR_ALL'], toggle = ['HYDR_USE', 'open', [ True ] ] )
    add_waters = self.createLine( [ 'widget', 'ADD_WATERS', 'label', 'Add waters' ] )
    self.createLine( [ 'label', '&nbsp;and then perform further ', 'widget', 'NCYCLES_AFTER_ADD_WATERS', 'label', ' refinement cycles' ], toggle = [ 'ADD_WATERS','open', [ True ] ], appendLine=add_waters  )
    self.createLine( [ 'widget', 'USE_TWIN', 'label', 'Crystal is twinned' ] )
    self.createLine( [ 'label', '<i>Warning: Intensities should be given for twin refinement. Using amplitudes is suboptimal.</i>' ], toggleFunction=[self.ToggleTwinSuboptimal, ['HKLIN_IS_I_SIGI', 'F_SIGF_OR_I_SIGI', 'HKLIN']])

    """self.createLine( [ 'label', '' ], toggleFunction=[self.ToggleTwinNotAvailable,['USEANOMALOUSFOR','HKLIN']])
    msg11 = 'Twin refinement not available for anomalous data.'
    line = self.createLine( [ 'label', msg11 ], toggleFunction=[self.ToggleTwinNotAvailable,['USEANOMALOUSFOR','HKLIN']])
    sizeFix = QtWidgets.QSizePolicy.Fixed
    btn = QtWidgets.QPushButton('help', line)
    btn.setSizePolicy(sizeFix, sizeFix)
    btn.released.connect(self.twinHelpPressed)
    line.layout().insertWidget(1, btn)"""

    self.closeSubFrame()
    
    #folder = self.openFolder(folderFunction='controlParameters',title='Options', drawFolder=self.drawControlParameters)
    #folder = self.openFolder(folderFunction='controlParameters',title='Advanced Options', drawFolder=self.drawAdvancedOptions)
    folder = self.openFolder(folderFunction='controlParameters',title='Parameterisation', drawFolder=self.drawParameters)
    folder = self.openFolder(folderFunction='controlParameters',title='Restraints', drawFolder=self.drawRestraints)
    #folder = self.openFolder(folderFunction='controlParameters',title='Output', drawFolder=self.drawOutput) # MM
    folder = self.openFolder(folderFunction='controlParameters',title='Advanced')
    self.drawAdvanced() # small change introduced to allow for automatically loading a keyword file in the 'advanced' tab
    self.setProsmartProteinMode()
    self.setProsmartNucleicAcidMode()
    self.setLibgMode()
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
    self.createLine( [ 'widget', 'B_REFINEMENT_MODE', 'label','ADPs'] )
    self.closeSubFrame()
    #self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']] )
    #self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    #self.closeSubFrame()
    
    self.createLine( [ 'subtitle', 'Scaling'] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'NO_SOLVENT', 'label', 'Do not consider bulk solvent contribution' ])
    #self.createLine( [ 'widget', 'SCALE_TYPE', 'label', 'solvent scaling, with', 'widget', 'SOLVENT_MASK_TYPE', 'label', 'solvent mask' ] )
    #self.createLine( [ 'widget', 'SCALE_TYPE', 'label', 'solvent scaling, with', 'widget', 'SOLVENT_MASK_TYPE', 'label', 'solvent mask' ], toggleFunction=[self.CheckScaleType, ['SCALE_TYPE']] )
    #self.createLine( [ 'widget', 'SOLVENT_ADVANCED', 'label', 'Use custom solvent mask parameters' ], toggleFunction=[self.ToggleUseSolventMask, ['SOLVENT_MASK_TYPE']] )
    #self.createLine( [ 'label', indent, 'label', 'Increase VDW radius of non-ion atoms by', 'widget', 'SOLVENT_VDW_RADIUS' ], toggleFunction=[self.ToggleSolventAdvanced, ['SOLVENT_MASK_TYPE', 'SOLVENT_ADVANCED']] )
    #self.createLine( [ 'label', indent, 'label', 'Increase ionic radius of potential ion atoms by', 'widget', 'SOLVENT_IONIC_RADIUS' ], toggleFunction=[self.ToggleSolventAdvanced, ['SOLVENT_MASK_TYPE', 'SOLVENT_ADVANCED']] )
    #self.createLine( [ 'label', indent, 'label', 'Shrink the mask area by', 'widget', 'SOLVENT_SHRINK', 'label', 'after calculation' ], toggleFunction=[self.ToggleSolventAdvanced, ['SOLVENT_MASK_TYPE', 'SOLVENT_ADVANCED']] )
    self.closeSubFrame()

    """self.createLine( [ 'subtitle', 'Translation-Libration-Screw (TLS)'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']] )
    self.createLine( [ 'label', 'TLS parameters', 'widget', 'TLSMODE' ], toggleFunction=[self.ToggleTLSModeOff, ['TLSMODE']] )
    self.createLine( [ 'label', 'TLS parameters', 'widget', 'TLSMODE', 'stretch', 'label', 'Number of TLS refinement cycles:', 'widget', 'NTLSCYCLES' ], toggleFunction=[ self.ToggleTLSModeOn, ['TLSMODE']] )
    self.createLine( [ 'label', '(To create TLS groups, use the "Import and/or edit TLS definitions" task in the Refinement Task menu)'], toggleFunction=[self.ToggleTLSModeFile, ['TLSMODE']] )
    self.createLine( [ 'widget', '-browseDb', True, 'TLSIN'], toggleFunction=[self.ToggleTLSModeFile, ['TLSMODE']] )
    reset_bfac_tls = self.createLine( [ 'widget', 'BFACSETUSE', 'label', 'Reset all B-factors at start' ], toggleFunction=[self.ToggleTLSModeOn, ['TLSMODE']])
    self.createLine( [ 'label', '&nbsp;to fixed value:', 'widget', 'BFACSET' ], toggle = ['BFACSETUSE', 'open', [ True ] ], appendLine=reset_bfac_tls )
    self.createLine( [ 'widget', 'TLSOUT_ADDU', 'label', 'Add TLS contribution to output B-factors (only for analysis and deposition)' ], toggleFunction=[self.ToggleTLSModeOn, ['TLSMODE']] )
    self.container.inputData.TLSIN.dataChanged.connect(self.checkAllowUndefined)
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']] )
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']] )
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()

    #self.createLine( [ 'label', 'TLS input file has been provided. Number of TLS refinement cycles:', 'stretch', 'widget', 'NTLSCYCLES' ], toggleFunction=[self.ToggleTLS, ['TLSIN']])
    #tls_cycles = self.createLine( [ 'widget', 'AUTOTLS', 'label', 'Use TLS parameters' ], toggleFunction=[self.ToggleTLSNot, ['TLSIN']])
    #self.createLine( [ 'label', 'Number of TLS refinement cycles:', 'widget', 'NTLSCYCLES_AUTO' ], toggleFunction=[self.ToggleTLSUsed, ['TLSIN','AUTOTLS']], appendLine=tls_cycles)
    #reset_bfac_tls = self.createLine( [ 'label', indent, 'widget', 'TLSBFACSETUSE', 'label', 'Reset all B-factors at start' ], toggleFunction=[self.ToggleTLSUsed, ['TLSIN','AUTOTLS']])
    #self.createLine( [ 'label', '&nbsp;to fixed value:', 'widget', 'TLSBFACSET' ], toggle = ['TLSBFACSETUSE', 'open', [ True ] ], appendLine=reset_bfac_tls )"""
    
    """self.createLine( [ 'subtitle', 'Rigid body groups'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']] )
    self.createLine( [ 'widget', 'RIGID_BODY_SELECTION' ] )
    self.closeSubFrame()
    self.getWidget('RIGID_BODY_SELECTION').setMinimumHeight(200)
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']] )
    self.createLine( [ 'label', '<i>Only available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()"""
    
    self.createLine( [ 'subtitle', 'Conformer groups and occupancy refinement'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']] )
    self.createLine( [ 'widget', 'OCCUPANCY_GROUPS', 'label', 'Specify partial occupancy groups (alternative conformers)' ] )
    self.createLine( [ 'widget', 'OCCUPANCY_SELECTION' ], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'OCCUPANCY_COMPLETE', 'label', 'Specify overlapping alternative conformer groups (constrain occupancies to sum to one)' ], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'OCCUPANCY_COMPLETE_TABLE' ], toggleFunction=[self.ToggleOccComplete, ['OCCUPANCY_GROUPS','OCCUPANCY_COMPLETE']] )
    self.createLine( [ 'widget', 'OCCUPANCY_INCOMPLETE', 'label', 'Specify overlapping alternative conformer groups (occupancies sum to less than one)' ], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
    self.createLine( [ 'widget', 'OCCUPANCY_INCOMPLETE_TABLE' ], toggleFunction=[self.ToggleOccIncomplete, ['OCCUPANCY_GROUPS','OCCUPANCY_INCOMPLETE']] )
    self.createLine( [ 'widget', 'OCCUPANCY_REFINEMENT', 'label', 'Perform refinement of atomic occupancies' ], toggle = ['OCCUPANCY_GROUPS', 'open', [ True ] ] )
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
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
    
    self.createLine( [ 'subtitle', 'Weights'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    auto_weight = self.createLine( [ 'label', 'Weight restraints versus experimental data using', 'widget', 'WEIGHT_OPT', 'label', 'weight'] )
    self.createLine( [ 'label', ':', 'widget', 'WEIGHT' ], toggleFunction=[self.ToggleWeightAuto, ['WEIGHT_OPT']], appendLine=auto_weight )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Non-Crystallographic Symmetry (NCS)'] )
    #self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOn,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    self.createLine( [ 'widget', 'USE_NCS', 'label', 'Use local non-crystallographic symmetry (NCS) restraints' ] )
    #self.createLine( [ 'label', indent+'Use automatic', 'widget', 'NCS_TYPE', 'label', 'NCS restraints' ], toggle = ['USE_NCS', 'open', [ True ] ] )
    self.closeSubFrame()
    #self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    #self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    #self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOff,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Covalent links'] )
    #self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'FIND_LINKS', 'label', 'Detect and apply covalent linkages based on the current atomic coordinates' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Jelly-body'] )
    #self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleJellyOn,['UNRESTRAINED', 'FIX_XYZ']])
    use_jellybody = self.createLine( [ 'widget', 'USE_JELLY', 'label', 'Use jelly-body restraints' ] )
    self.createLine( [ 'label', '&nbsp;with sigma:', 'widget', 'JELLY_SIGMA', 'label', 'and max distance:', 'widget', 'JELLY_DIST' ], toggle = ['USE_JELLY', 'open', [ True ] ], appendLine=use_jellybody)
    self.createLine( [ 'widget', 'JELLY_ONLY', 'label', 'Jelly body refinement only' ], toggle = ['USE_JELLY', 'open', [ True ] ] )
    self.closeSubFrame()
    #self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    #self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    #self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleJellyOff,['UNRESTRAINED', 'FIX_XYZ']])
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'MetalCoord External Restraints for Metals'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOn,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    self.createLine( [ 'widget', 'metalCoordPipeline.RUN_METALCOORD', 'label', 'Apply MetalCoord restraints for metal sites:' ] , toggle = ['metalCoordPipeline.RUN_METALCOORD', 'open', [ False ] ])
    self.createLine( [ 'widget', 'metalCoordPipeline.RUN_METALCOORD', 'label', 'Apply MetalCoord restraints for metal sites:' , 'stretch' , 'widget', 'metalCoordPipeline.GENERATE_OR_USE'] , toggle = ['metalCoordPipeline.RUN_METALCOORD', 'open', [ True ] ])
    self.createLine( [ 'widget', 'metalCoordPipeline.METALCOORD_RESTRAINTS'] , toggleFunction=[self.ToggleMetalCoordUse, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE']])
    if self.widget.subFrame is not None:
        self.currentLayout = self.widget.subFrame.layout()
    else:
        self.currentLayout = self.widget.currentFolderLayout
    self.ligands_checkboxes = ChoiceButtons()
    self.currentLayout.addWidget(self.ligands_checkboxes)
    self.updateMonomersWithMetalsWidget()
    self.ligands_checkboxes.clickedSignal.connect(self.updateMonomersWithMetalsSelection)
    self.createLine( [ 'widget', 'metalCoordPipeline.TOGGLE_ADVANCED', 'label', 'Show advanced options' ], toggleFunction=[self.ToggleMetalCoordGenerate, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE']] )
    self.createLine( [ 'label', 'Link records to metal sites in the atomic model:', 'stretch', 'widget', 'metalCoordPipeline.LINKS' ], toggle = ['metalCoordPipeline.TOGGLE_ADVANCED', 'open', [ True ] ] )
    #self.createLine( [ 'widget', 'metalCoordWrapper.KEEP_LINKS', 'label', 'Do not delete the link records to metal sites which are already present in the atomic model' ], toggle = ['metalCoordPipeline.UPDATE_LINKS', 'open', [ True ] ] )
    self.createLine( [ 'label', 'Distance threshold: (range 0-1)<br/><i>A threshold d to select atoms is (r<sub>1</sub> + r<sub>2</sub>)*(1 + d) where r<sub>1</sub> and r<sub>2</sub> are covalent radii.</i>', 'stretch', 'widget', 'metalCoordWrapper.DISTANCE_THRESHOLD' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    self.createLine( [ 'label', 'Maximum coordination number:', 'stretch', 'widget', 'metalCoordWrapper.MAXIMUM_COORDINATION_NUMBER' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    self.createLine( [ 'label', 'Procrustes distance threshold: (range 0-1)', 'stretch', 'widget', 'metalCoordWrapper.PROCRUSTES_DISTANCE_THRESHOLD' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    self.createLine( [ 'label', 'Minimum sample size for statistics:', 'stretch', 'widget', 'metalCoordWrapper.MINIMUM_SAMPLE_SIZE' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    self.createLine( [ 'widget', 'metalCoordWrapper.USE_PDB', 'label', 'Use COD structures based on the input PDB/mmCIF coordinates' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    self.createLine( [ 'widget', 'metalCoordWrapper.IDEAL_ANGLES', 'label', 'Provide only ideal bond angles' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    self.createLine( [ 'widget', 'metalCoordWrapper.SIMPLE', 'label', 'Simple distance based filtering' ], toggleFunction=[self.ToggleMetalCoordGenerateAdvanced, ['metalCoordPipeline.RUN_METALCOORD', 'metalCoordPipeline.GENERATE_OR_USE', 'metalCoordPipeline.TOGGLE_ADVANCED']] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRestraintsOff,['UNRESTRAINED', 'FIX_XYZ', 'JELLY_ONLY']])
    self.createLine( [ 'label', '<i>Not available.</i>' ] )
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
    #self.createLine( [ 'widget', '-browseDb', True, 'prosmartProtein.REFERENCE_MODEL' ] )
    self.createLine( [ 'widget', '-browseDb', True, 'prosmartProtein.REFERENCE_MODELS' ] )
    self.createLine( [ 'label', 'Use', 'widget', 'prosmartProtein.ALL_BEST', 'label', 'chain(s) from the reference model(s).', 'stretch', 'label', 'Minimum sequence identity:', 'widget', 'prosmartProtein.SEQID', 'label', '%' ] )
    self.createLine( [ 'label', 'Generate restraints between', 'widget', 'prosmartProtein.SIDE_MAIN', 'label', 'atom-pairs.', 'stretch', 'label', 'Interatomic distance range:', 'widget', 'prosmartProtein.RMIN', 'label', 'to', 'widget', 'prosmartProtein.RMAX', 'label', 'angstroem' ] )
    ###self.createLine( [ 'label', 'Apply restraints up to distance', 'widget', 'prosmartProtein.DMAX', 'label', 'angstroem with weight', 'widget', 'prosmartProtein.WEIGHT', 'label', 'and robustness parameter (alpha)', 'widget', 'prosmartProtein.ALPHA', 'stretch', 'label', 'Show advanced options', 'widget', 'prosmartProtein.ADVANCED' ])
    self.createLine( [ 'label', 'Apply restraints up to distance', 'widget', 'prosmartProtein.DMAX', 'label', 'angstroem with robustness parameter (alpha)', 'widget', 'prosmartProtein.ALPHA' , 'stretch', 'label', 'Show advanced options', 'widget', 'prosmartProtein.ADVANCED' ])
    ##self.createLine( [ 'label', 'Apply restraints with weight', 'widget', 'prosmartProtein.WEIGHT', 'label', 'and robustness parameter (alpha)', 'widget', 'prosmartProtein.ALPHA', 'stretch' , 'label', 'Show advanced options', 'widget', 'prosmartProtein.ADVANCED' ])
    self.createLine( [ 'label', 'Minimum and maximum sigma', 'widget', 'prosmartProtein.SGMN', 'widget', 'prosmartProtein.SGMX'], toggle = ['prosmartProtein.ADVANCED', 'open', [ True ] ] )    
    self.createLine( [ 'widget', 'prosmartProtein.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high ADPs.' ] , toggleFunction=[self.hideProsmartProteinBfac, ['prosmartProtein.ADVANCED','prosmartProtein.TOGGLE_BFAC']])
    ##self.createLine( [ 'widget', 'prosmartProtein.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high B-factors.' ] , toggleFunction=[self.hideProsmartProteinBfac, 'prosmartProtein.TOGGLE_BFAC'])
    self.createLine( [ 'widget', 'prosmartProtein.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high ADPs.', 'stretch', 'label', 'Maximum B-factor: median plus ', 'widget', 'prosmartProtein.BFAC', 'label', 'x interquartile range' ] , toggleFunction=[self.showProsmartProteinBfac, ['prosmartProtein.ADVANCED','prosmartProtein.TOGGLE_BFAC']])
    ##self.createLine( [ 'widget', 'prosmartProtein.TOGGLE_BFAC', 'label', 'Remove restraints where homologue has high B-factors.', 'stretch', 'label', 'Maximum B-factor: median plus ', 'widget', 'prosmartProtein.BFAC', 'label', 'x interquartile range' ] , toggleFunction=[self.showProsmartProteinBfac, ['prosmartProtein.TOGGLE_BFAC']])
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

    '''
    self.createLine( [ 'subtitle', 'LibG External Restraints for Nucleic Acids'] )
    if self.isEditable():
       self.container.libg.TOGGLE.dataChanged.connect( self.setLibgMode)
       self.container.controlParameters.REFINEMENT_MODE.dataChanged.connect(self.setLibgMode)
    self.openSubFrame(frame=[True], toggle = ['libg.MODE', 'open', [ 'DISABLED' ] ] )
    self.createLine( [ 'label', '<i>Specify atomic model before setting up external restraints</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['libg.MODE', 'open', [ 'NONUCLEICACID' ] ] )
    self.createLine( [ 'label', '<i>Input atomic model contains no nucleotide chains</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['prosmartProtein.MODE', 'open', [ 'RIGIDMODE' ] ] )
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['libg.MODE', 'open', [ 'UNSELECTED' ] ] )
    self.createLine( [ 'widget', 'libg.TOGGLE', 'label', 'Generate restraints for nucleic acids' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggle = ['libg.MODE', 'open', [ 'SELECTED' ] ] )
    self.createLine( [ 'widget', 'libg.TOGGLE', 'label', 'Generate restraints for nucleic acids:', 'widget', 'libg.OPTION' ] )
    self.createLine( [ 'label', indent, 'widget', 'libg.BP', 'label', 'Base-pairs' ], toggleFunction=[self.showLibgOptions, ['libg.OPTION']])
    self.createLine( [ 'widget', 'libg.ADVANCED', 'label', 'Show advanced options' ])
    self.createLine( [ 'label', 'Additional LibG keywords:', 'widget', 'libg.KEYWORDS' ], toggle = ['libg.ADVANCED', 'open', [ True ] ] )
    self.closeSubFrame()
    '''

    """self.createLine( [ 'subtitle', 'Platonyzer Metal Site Restraints'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    use_platonyzer = self.createLine( [ 'widget', 'platonyzer.TOGGLE', 'label', 'Use Platonyzer restraints' ] )
    self.createLine( [ 'label', '&nbsp;for', 'widget', 'platonyzer.MODE', 'label', 'ions' ], toggle = ['platonyzer.TOGGLE', 'open', [ True ] ], appendLine=use_platonyzer)
    self.createLine( [ 'widget', 'platonyzer.RM_VDW', 'label', 'Remove VDW restraints for octohedral ions' ], toggleFunction=[self.TogglePlatonyzerNAMG,['platonyzer.TOGGLE','platonyzer.MODE']] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()"""

    """self.createLine( [ 'subtitle', 'Additional Restraints'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    self.createLine( [ 'widget', '-browseDb', True, 'EXTERNAL_RESTRAINTS_FILE' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()"""
    
    """self.createLine( [ 'subtitle', 'Advanced Options'] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    self.createLine( [ 'widget', 'MAKE_LINK', 'label', 'Detect and apply covalent linkages based on the current atomic coordinates' ] )
    self.createLine( [ 'widget', 'OVERRIDE_LINK', 'label', 'Ignore all LINK records in the input model' ], toggle = ['MAKE_LINK', 'open', [ True ] ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()"""

    self.createLine( [ 'subtitle', 'ADP Restraints'] )
    self.openSubFrame()
    self.createLine( [ 'label', 'ADP restraint weight:', 'stretch', 'widget', 'ADPR_WEIGHT' ] )
    self.createLine( [ 'label', 'Maximum distance for ADP restraint:', 'stretch', 'widget', 'MAX_DIST_FOR_ADP_RESTRAINT' ] )
    self.createLine( [ 'label', 'ADP restraint power:', 'widget', 'stretch', 'ADP_RESTRAINT_POWER' ] )
    self.createLine( [ 'widget', 'ADP_RESTRAINT_NO_LONG_RANGE', 'label', 'No long range for ADP restraint' ] )
    #self.createLine( [ 'label', 'ADP restraint mode:', 'stretch', 'widget', 'ADP_RESTRAINT_MODE' ] )
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Infrequently Used Options'] )
    self.openSubFrame()
    self.createLine( [ 'widget', 'UNRESTRAINED', 'label', 'No positional restraints' ] )
    self.createLine( [ 'widget', 'FIX_XYZ', 'label', 'Fix coordinates' ] )
    self.closeSubFrame()
    return

  def showProsmartProteinBfac(self):
     # if self.container.prosmartProtein.ADVANCED and self.container.prosmartProtein.TOGGLE_BFAC:
     if self.container.prosmartProtein.TOGGLE_BFAC:
         return True
     return False

  def showProsmartNucleicAcidBfac(self):
     # if self.container.prosmartNucleicAcid.ADVANCED and self.container.prosmartNucleicAcid.TOGGLE_BFAC:
     if self.container.prosmartNucleicAcid.TOGGLE_BFAC:
         return True
     return False
         
  def hideProsmartProteinBfac(self):
     if self.container.prosmartProtein.ADVANCED and not self.container.prosmartProtein.TOGGLE_BFAC:
        ##if not self.container.prosmartProtein.TOGGLE_BFAC:
        return True
     return False

  def hideProsmartNucleicAcidBfac(self):
     if self.container.prosmartNucleicAcid.ADVANCED and not self.container.prosmartNucleicAcid.TOGGLE_BFAC:
        ##if not self.container.prosmartNucleicAcid.TOGGLE_BFAC:
        return True
     return False
  
  def showLibgOptions(self):
     if str(self.container.libg.OPTION) == 'MANUAL':
        return True
     return False

  def drawOutput( self ):
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

    """self.createLine( [ 'subtitle', 'Output Options' ], toggleFunction=[self.ToggleNeutronModeOff,['SCATTERING_FACTORS']] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleNeutronModeOff,['SCATTERING_FACTORS']])
    self.createLine( [ 'label', 'Output calculated riding hydrogens to file', 'widget', 'OUTPUT_HYDROGENS' ] )
    self.closeSubFrame()"""

    """self.createLine( [ 'subtitle', 'Map Calculation' ] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'MAP_SHARP', 'label', 'Perform map sharpening when calculating maps' ] )
    use_mapsharp = self.createLine( [ 'label', indent, 'widget', 'MAP_SHARP_CUSTOM', 'label', 'Use custom sharpening parameter (B-factor)' ], toggle = ['MAP_SHARP', 'open', [ True ] ] )
    self.createLine( [ 'label', ':', 'widget', 'BSHARP' ], toggle = ['MAP_SHARP_CUSTOM', 'open', [ True ] ], appendLine=use_mapsharp)
    self.createLine( [ 'label', '' ] )

    self.createLine( [ 'label', '<i>Anomalous maps will not be generated - anomalous data not provided in Input Data</i>' ] , toggle = ['USEANOMALOUS', 'open', [ False ] ] )
    self.createLine( [ 'label', '<i>Anomalous maps will be generated</i>' ] , toggle = ['USEANOMALOUS', 'open', [ True ] ] )
    self.closeSubFrame()"""

    """self.createLine( [ 'subtitle', 'Validation and Analysis' ] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOff,['REFINEMENT_MODE']])
    self.createLine( [ 'widget', 'VALIDATE_IRIS', 'label', 'Generate Iris report' ] )
    self.createLine( [ 'widget', 'VALIDATE_BAVERAGE', 'label', 'Analyse B-factor distributions' ] )
    self.createLine( [ 'widget', 'VALIDATE_RAMACHANDRAN', 'label', 'Generate Ramachandran plots' ] )
    self.createLine( [ 'widget', 'VALIDATE_MOLPROBITY', 'label', 'Run MolProbity to analyse geometry' ] )
    #self.createLine( [ 'widget', 'RUN_MOLPROBITY', 'label', 'Run standalone MolProbity (to be deprecated)' ] )
    self.closeSubFrame()
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleRigidModeOn,['REFINEMENT_MODE']])
    self.createLine( [ 'label', '<i>Not available in Rigid Body mode.</i>' ] )
    self.closeSubFrame()"""
    return

  def drawAdvanced( self ):
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
    #- REMOVED: self.createLine( [ 'label','Exit if new ligand encountered','stretch','widget', 'MAKE_NEW_LIGAND_EXIT'] )

    """self.createLine( [ 'label', indent+'<i>(see Parameterisation tab for H/D fraction refinement)</i>' ], toggleFunction=[self.ToggleNeutronModeOn, ['SCATTERING_FACTORS']], appendLine=scattering_factors )
    self.createLine( [ 'subtitle', 'Neutron refinement options'], toggleFunction=[self.ToggleNeutronModeOn,['SCATTERING_FACTORS']] )
    self.openSubFrame(frame=[True], toggleFunction=[self.ToggleNeutronModeOn,['SCATTERING_FACTORS']] )
    self.createLine( [ 'widget', 'HYDR_USE', 'label', 'Use hydrogens during refinement'], toggle = ['HYDR_USE', 'open', [ False ] ])
    self.createLine( [ 'widget', 'HYDR_USE', 'label', 'Use hydrogens during refinement', 'widget', 'HYDR_ALL'], toggle = ['HYDR_USE', 'open', [ True ] ] )
    use_h = self.createLine( [ 'widget', 'H_REFINE', 'label', 'Refine hydrogen positions' ], toggle = ['HYDR_USE', 'open', [ True ] ] )
    self.createLine( [ 'label', 'for', 'widget', 'H_REFINE_SELECT' ], toggleFunction=[self.ToggleH_REFINE,['HYDR_USE','H_REFINE']], appendLine=use_h)
    self.createLine( [ 'widget', 'H_TORSION', 'label', 'Use hydrogen torsion angle restraints' ], toggle = ['HYDR_USE', 'open', [ True ] ] )
    use_hd = self.createLine( [ 'widget', 'HD_FRACTION', 'label', 'Refine hydrogen/deuterium fractions' ], toggle = ['HYDR_USE', 'open', [ True ] ] )
    self.createLine( [ 'label', 'for', 'widget', 'HD_FRACTION_TYPE' ], toggleFunction=[self.ToggleNeutronModeHD,['HYDR_USE','HD_FRACTION']], appendLine=use_hd)
    self.createLine( [ 'label', indent, 'label', 'Initialise H/D fractions', 'widget', 'HD_INIT' ], toggleFunction=[self.ToggleNeutronModeHD_YES,['HYDR_USE','HYDR_ALL','HD_FRACTION']])
    self.createLine( [ 'label', indent, 'label', 'Initialise H/D fractions', 'widget', 'HD_INIT_HALL' ], toggleFunction=[self.ToggleNeutronModeHD_ALL,['HYDR_USE','HYDR_ALL','HD_FRACTION']])
    self.closeSubFrame()"""

    custom_res = self.createLine( [ 'widget', 'RES_CUSTOM', 'label', 'Use custom resolution limits' ] )
    self.createLine( [ 'label', indent+indent+'highest (d<sub>min</sub>):', 'widget', 'RES_MIN', 'label', ' lowest (d<sub>max</sub>):', 'widget', 'RES_MAX' ], toggle = ['RES_CUSTOM', 'open', [ True ] ], appendLine=custom_res )
    self.createLine( [ 'label', 'FreeR flag number for test set:' , 'widget', 'FREERFLAG_NUMBER'])
    scattering_factors = self.createLine( [ 'label', 'Diffraction experiment type:', 'widget', 'SCATTERING_FACTORS' ] )
    #if self.isEditable():
    #    self.container.controlParameters.SCATTERING_FACTORS.dataChanged.connect( self.ExperimentChanged)
    #self.createLine( [ 'label', indent+indent+'Form factor calculation:', 'widget', 'SCATTERING_ELECTRON' ], toggleFunction=[self.ToggleElectronDiffraction, ['SCATTERING_FACTORS']], appendLine=scattering_factors )
    self.createLine( [ 'widget', 'H_REFINE', 'label', 'Refine hydrogen positions' ], toggle = ['HYDR_USE', 'open', [ True ] ] )

    self.createLine( [ 'subtitle', 'Structure model modification before refinement' ] )
    self.openSubFrame(frame=[True])
    reset_bfac = self.createLine( [ 'widget', 'BFACSETUSE', 'label', 'Reset all ADPs at start' ])
    self.createLine( [ 'label', '&nbsp;to fixed value:', 'widget', 'BFACSET' ], toggle = ['BFACSETUSE', 'open', [ True ] ], appendLine=reset_bfac )
    randomize = self.createLine( [ 'widget', 'RANDOMIZEUSE', 'label', 'Shake coordinates at start' ])
    self.createLine( [ 'label', '&nbsp;with specified RMSD:', 'widget', 'RANDOMIZE' ], toggle = ['RANDOMIZEUSE', 'open', [ True ] ], appendLine=randomize )
    self.closeSubFrame()

    self.createLine( [ 'widget', 'USE_WORK_IN_EST', 'label', 'Use work reflections in ML parameter estimates' ] )
    self.createLine( [ 'widget', 'KEEP_CHARGES', 'label', 'Keep charges, i.e. use scattering factor for charged atoms where relevant' ] )
    #self.createLine( [ 'widget', 'REFMAC_CLEANUP', 'label', 'Clean up intermediate files at end of job' ] )

    self.createLine( [ 'subtitle', 'Additional keywords'] )
    # self.createLine( [ 'label', '<i>Keywords specified below will overwrite options which were set elsewhere.</i>'] )
    # self.createLine( [ 'widget', '-guiMode','multiLine','EXTRAREFMACKEYWORDS' ] )
    self.createLine( [ 'widget', '-browseDb', True, 'SERVALCAT_KEYWORD_FILE' ] )
    self.createLine( [ 'label', 'Extra servalcat command line options:', 'widget', 'EXTRA_SERVALCAT_OPTIONS' ] )
    self.getWidget('EXTRA_SERVALCAT_OPTIONS').setFixedWidth(400)

    self.createLine( [ 'subtitle', 'Validation and Analysis' ] )
    self.openSubFrame(frame=[True])
    self.createLine( [ 'widget', 'VALIDATE_IRIS', 'label', 'Generate Iris validation report' ] )
    # self.createLine( [ 'widget', 'VALIDATE_BAVERAGE', 'label', 'Analyse B-factor distributions' ] )
    self.createLine( [ 'widget', 'VALIDATE_RAMACHANDRAN', 'label', 'Generate Ramachandran plots' ] )
    self.createLine( [ 'widget', 'VALIDATE_MOLPROBITY', 'label', 'Run MolProbity to analyse geometry' ] )
    # self.closeSubFrame()

    # self.createLine( [ 'subtitle', 'Monitoring' ] )
    # self.openSubFrame(frame=[True])
    # TO DO - on/off
    self.createLine( [ 'label', 'Minimum deviation of atom coordinates to be reported:', 'stretch', 'widget', 'monitor.MIN_COORDDEV' ] )
    self.createLine( [ 'label', 'Minimum deviation of B-values to be reported:', 'stretch', 'widget', 'monitor.MIN_ADPDEV' ] )
    self.createLine( [ 'label', 'Atoms with a B-value lower than <i>the first quartile - factor * interquartile_range</i><br />or higher than <i>the third quartile + factor * interquartile_range</i> to be reported. Factor:', 'stretch', 'widget', 'monitor.ADP_IQR_FACTOR' ] )
    self.closeSubFrame()
    return

#-  --------------------          --------------------          --------------------

  def anomalousAvailable(self):
    if not self.container.inputData.HKLIN.isSet(): return False
    if not self.isEditable():
       if self.container.controlParameters.USEANOMALOUS:
          return True #only display after job is run if USEANOMALOUS has already been set to True
       else:
          return False
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

  #def anomalousMapUsed(self):
  #  if not self.container.inputData.HKLIN.isSet(): return False
  #  if not self.container.controlParameters.USEANOMALOUSFOR.isSet(): return False
  #  if str(self.container.controlParameters.USEANOMALOUSFOR) == 'NOTHING': return False
  #  return True

#def anomalousMapIgnored(self):
#   if self.anomalousAvailable():
#      if self.container.controlParameters.USEANOMALOUSFOR.isSet():
#         if str(self.container.controlParameters.USEANOMALOUSFOR) == 'NOTHING': return True
#   return False


  @QtCore.Slot()
  def hklinChanged(self):
    self.getWavelength()
    self.getObsType()
  #self.setTwinMode()

  #def setTwinMode(self):
  #  if self.container.inputData.HKLIN.isSet():
  #     columnLabelsInFile = [column.columnLabel.__str__() for column in self.container.inputData.HKLIN.fileContent.listOfColumns]
  #     if not 'I' in columnLabelsInFile and not 'Iplus' in columnLabelsInFile:
  #         self.container.controlParameters.TWIN_TYPE.setQualifiers({'enumerators':['F'],'menuText':['SF Amplitudes (F)']})
  #         self.container.controlParameters.TWIN_TYPE.set('F')
  #     else:
  #         self.container.controlParameters.TWIN_TYPE.setQualifiers({'enumerators':['I','F'],'menuText':['Intensities (I)','SF Amplitudes (F)']})
  #         self.container.controlParameters.TWIN_TYPE.set('I')
  #     try:
  #         self.getWidget('TWIN_TYPE').populateComboBox(self.container.controlParameters.TWIN_TYPE)
  #         self.getWidget('TWIN_TYPE').updateViewFromModel()
  #     except:
  #         pass
  #     self.validate()
  
  @QtCore.Slot()
  def ExperimentChanged(self):
    if str(self.container.controlParameters.SCATTERING_FACTORS) == "NEUTRON":
      self.container.controlParameters.H_REFINE.set(True)
    else:
      self.container.controlParameters.H_REFINE.set(False)
    return
    
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
     self.setLibgMode()
     self.getMonomersWithMetals()
     return

  @QtCore.Slot()
  def setProsmartProteinMode(self):
     self.container.prosmartProtein.CHAINLIST_1.setQualifiers({'allowUndefined':True})
     #self.container.prosmartProtein.REFERENCE_MODEL.setQualifiers({'allowUndefined':True})
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
                 #self.container.prosmartProtein.REFERENCE_MODEL.setQualifiers({'allowUndefined':False})
                 self.container.prosmartProtein.REFERENCE_MODELS.setQualifiers({'listMinLength':1})
              else:
                 self.container.prosmartProtein.MODE.set('UNSELECTED')
              self.validate()
              return
        self.container.prosmartProtein.MODE.set('NOPROTEIN')
     else:
        self.container.prosmartProtein.MODE.set('DISABLED')
     self.validate()
     return

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
     return

  @QtCore.Slot()
  def setLibgMode(self):
     if self.ToggleRigidModeOn():
        self.container.libg.MODE.set('RIGIDMODE')
     elif self.container.inputData.XYZIN.isSet():
        if self.container.inputData.NUCLEOTIDE_CHAINS.isSet():
           if self.container.inputData.NUCLEOTIDE_CHAINS:
              if self.container.libg.TOGGLE:
                 self.container.libg.MODE.set('SELECTED')
              else:
                 self.container.libg.MODE.set('UNSELECTED')
              self.validate()
              return
        self.container.libg.MODE.set('NONUCLEICACID')
     else:
        self.container.libg.MODE.set('DISABLED')
     self.validate()
     return

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
         #print "Model content: "+str(chain_list)
     return chain_list

#      from iotbx.pdb import hierarchy
#      pdb_in = hierarchy.input(file_name=self.container.inputData.XYZIN.fullPath.__str__())
#      for chain in pdb_in.hierarchy.only_model().chains():
#         if chain.is_protein():
#            print "protein: "+chain.id
#         elif chain.is_na():
#            print "dna/rna: "+chain.id
#         else:
#            print "neither: "+chain.id

  #def MergeDictionaries(self):
  #   print 'Dictionaries are not merged here - its now done when the job is executed...'

  def getWavelength(self):
    return  # MM
    #if self.container.inputData.HKLIN.isSet():
    #    wavelengths = self.container.inputData.HKLIN.fileContent.getListOfWavelengths()
    #    if len(wavelengths)>0: self.container.controlParameters.WAVELENGTH = round(wavelengths[-1],3)
    #    try:
    #        self.getWidget('WAVELENGTH').updateViewFromModel()
    #    except:
    #        print('prosmart_refmac - WAVELENGTH widget  updateViewFromModel failed')

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
      invalidElements = super(Cservalcat_xtal_pipe, self).isValid()
      #Check whether invocation is from runTask
      import traceback
      functionNames = [a[2] for a in traceback.extract_stack()]

      #if str(self.container.controlParameters.DATA_METHOD) == 'xtal':
      if True:
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
               #from PyQt4.QtCore import *
               msg = QMessageBox()
               msg.setIcon(QMessageBox.Question)

               msg.setText("Warning")
               msg.setInformativeText("You are trying to execute refinement without selecting a Free-R data item")
               msg.setWindowTitle("Free-R warning")
               msg.setDetailedText("While refinement without FreeR is reasonable under some circumstances, it is generally discouraged because it can provide a misleading impression of progress in refinement where over-fitting of the data can occur")
               msg.setStandardButtons(QMessageBox.Ignore | QMessageBox.Cancel)
               #msg.buttonClicked.connect(msgbtn)
               retval = msg.exec_()
               #print "value of pressed message box button:", retval, 0x00100000, 0x00400000
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


         return invalidElements

