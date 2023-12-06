"""
     tasks/chainsaw/CTaskChainsaw.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class CTaskChainsaw(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'chainsaw'
  TASKVERSION = 0.1
  TASKMODULE='bioinformatics'
  TASKTITLE='Truncate search model - CHAINSAW'
  SHORTTASKTITLE='CHAINSAW'
  DESCRIPTION='Truncate and renumber model prior to molecular replacement'
  MGDISPLAYFILES = ['XYZIN','XYZOUT']
  PERFORMANCECLASS = 'CAtomCountPerformance'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('chainsaw')
                        
    self.openFolder(folderFunction='inputData')
    self.createLine ( [ 'subtitle', 'Enter structure to be edited' ] )
    self.createLine ( [ 'tip','input search model from homologous protein',
                        'widget','XYZIN' ] )
    self.createLine ( [ 'subtitle', 'Enter sequence alignment of edited structure and MR target' ] )
    self.setToolTip('ALIGNIN','Sequence alignment of search model to target')
    self.createLine ( [ 'widget', 'ALIGNIN' ] )
    self.container.inputData.ALIGNIN.dataChanged.connect(self.alignmentChanged)
    self.createLine ( [ 'label','Identifier of the target sequence in this alignemnt','widget','-guiMode','combo','TARGETINDEX'] )
    self.alignmentChanged()
    
    self.openFolder(title='Simple Options')

    self.createLine( [ 'label', 'How severe is side-chain truncation for non-conserved residues: ' ] )
    self.setMenuText('MODE',{ 'MIXA':'truncate to CB',
                              'MIXS':'truncate to CG',
                              'MAXI':'truncate to last common atom'})
    self.createLine( [ 'widget', 'MODE',
                       ] )
  @QtCore.Slot()  
  def alignmentChanged(self):
      # Keep in sync with similar method in CTaskSculptor
      import os
      idList = []
      enumerators = []
      if self.container.inputData.ALIGNIN.exists():
        err,idList = self.container.inputData.ALIGNIN.bioGetSeqIdentifiers()
        for ii in range(len(idList)): enumerators.append(ii)
      self.container.controlParameters.TARGETINDEX.setQualifiers({'onlyEnumerators': True, 'enumerators':enumerators,'menuText':idList})
      self.getWidget('TARGETINDEX').populateComboBox(self.container.controlParameters.TARGETINDEX)
      self.container.controlParameters.TARGETINDEX.set(0)
      self.updateViewFromModel()
