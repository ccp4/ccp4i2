"""
     tasks/sculptor/CTaskSculptor.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.baselayer import QtCore

#-------------------------------------------------------------------
class CTaskSculptor(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'sculptor'
  TASKVERSION = 0.1
  TASKMODULE='bioinformatics'
  TASKTITLE='Truncate search model - SCULPTOR'
  SHORTTASKTITLE='SCULPTOR'
  MGDISPLAYFILES = ['XYZIN','XYZOUT']
  DESCRIPTION='Truncate model prior to molecular replacement'
  PERFORMANCECLASS = 'CAtomCountPerformance'
  ERROR_CODES = { 201 : { 'description' : 'ValueError parsing fasta format...try pushing alignment through clustalw to correct mismatched sequence lengths' },202 : { 'description' : 'Unknown problem with parsing fasta format...try pushing alignment through clustalw to correct mismatched sequence lengths' },}

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('sculptor')
    
                        
    self.openFolder(folderFunction='inputData')

    self.createLine ( [ 'subtitle', 'Template for modeling'])
    self.openSubFrame(frame=True)
    self.createLine ( [ 'advice','Template PDB to be used in modelling'] )
    self.createLine ( [ 'widget','XYZIN' ] )
    self.container.inputData.XYZIN.dataChanged.connect(self.xyzinChanged)
    self.createLine ( [ 'label','Chain ID in template to manipulate','widget', 'CHAINIDS' ],toggle=['ALIGNMENTORSEQUENCEIN','open',['SEQUENCE']] )
    self.xyzinChanged()
    self.closeSubFrame()
    
    self.createLine ( [ 'subtitle', 'Target sequence'])
    self.createLine ( ['label','Provide target sequence as :','widget','ALIGNMENTORSEQUENCEIN'])
    self.openSubFrame(frame=True, toggle=['ALIGNMENTORSEQUENCEIN','open',['ALIGNMENT']])
    self.createLine ( [ 'widget', 'ALIGNIN'] )
    self.container.inputData.ALIGNIN.dataChanged.connect(self.alignmentChanged)
    self.createLine ( [ 'label','Identifier of the target sequence in this alignemnt','widget','-guiMode','combo','TARGETINDEX'] )
    self.alignmentChanged()
    self.closeSubFrame()
    
    self.openSubFrame(frame=True, toggle=['ALIGNMENTORSEQUENCEIN','open',['SEQUENCE']])
    self.createLine ( [ 'widget', 'SEQUENCEIN'] )
    self.closeSubFrame()
    
    self.openFolder(title='Simple Options')

    self.setMenuText('PRUNING',{ 'schwarzenbacher':'mixed model of Schwarzenbacher et al.',
                              'similarity':'based on sequence similarity to target'})
    self.createLine( [ 'label', 'Pruning side-chains of search model: ',
                       'widget', 'PRUNING',
                       ] )

    self.setMenuText('BFACTOR',{ 'asa':'based on solvent accessibility',
                              'similarity':'based on sequence similarity to target',
                              'original':'keep the original'})
    self.createLine( [ 'label', 'Setting B factors of search model: ',
                       'widget', 'BFACTOR',
                       ] )

    self.container.inputData.ALIGNMENTORSEQUENCEIN.dataChanged.connect(self.alignmentOrSequenceChanged)
    self.alignmentOrSequenceChanged()

  @QtCore.Slot()
  def alignmentChanged(self):
      # Keep in sync with similar method in CTaskChainsaw
      import os
      idList = []
      enumerators = []
      if self.container.inputData.ALIGNIN.exists():
        err,idList = self.container.inputData.ALIGNIN.bioGetSeqIdentifiers()
        for ii in range(len(idList)): enumerators.append(ii)
      self.container.controlParameters.TARGETINDEX.setQualifiers({'onlyEnumerators': True, 'enumerators':enumerators,'menuText':idList})
      self.getWidget('TARGETINDEX').populateComboBox(self.container.controlParameters.TARGETINDEX)
      self.container.controlParameters.TARGETINDEX.set(0)
      self.getWidget('TARGETINDEX').updateViewFromModel()

  @QtCore.Slot()
  def xyzinChanged(self):
      import os
      if self.container.inputData.XYZIN.isSet() and os.path.isfile(self.container.inputData.XYZIN.__str__()):
        from core.CCP4ModelData import CPdbData
        aCPdbData = CPdbData()
        aCPdbData.loadFile(self.container.inputData.XYZIN.fullPath)
        enumerators = []
        menuText = []
        for chain in aCPdbData.composition.peptides:
            enumerators += [chain]
            menuText += [chain]
        if len(enumerators) > 0 and len(menuText) > 0:
            self.container.controlParameters.CHAINIDS.setQualifiers({'enumerators':enumerators})
            self.container.controlParameters.CHAINIDS.setQualifiers({'menuText':menuText})
            self.getWidget('CHAINIDS').populateComboBox(self.container.controlParameters.CHAINIDS)
            self.container.controlParameters.CHAINIDS.set(enumerators[0])
            self.updateViewFromModel()

  @QtCore.Slot()
  def alignmentOrSequenceChanged(self):
    if self.container.inputData.ALIGNMENTORSEQUENCEIN.__str__() == 'ALIGNMENT':
        self.container.inputData.ALIGNIN.setQualifiers({'allowUndefined':False})
        self.container.inputData.SEQUENCEIN.setQualifiers({'allowUndefined':True})
        self.container.inputData.SEQUENCEIN.unSet()
    else:
        self.container.inputData.ALIGNIN.setQualifiers({'allowUndefined':True})
        self.container.inputData.SEQUENCEIN.setQualifiers({'allowUndefined':False})
        self.container.inputData.ALIGNIN.unSet()

    self.getWidget('ALIGNIN').validate()
    self.getWidget('SEQUENCEIN').validate()

