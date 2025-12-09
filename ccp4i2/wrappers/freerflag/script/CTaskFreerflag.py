"""
     tasks/freerflag/CTaskFreerflag.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""
from ccp4i2.baselayer import QtCore
from qtgui.CCP4TaskWidget import CTaskWidget

class CTaskFreerflag(CTaskWidget):

  TASKNAME = 'freerflag'
  TASKVERSION = 0.0
  TASKMODULE= [ 'data_reduction', 'expt_data_utility' ]
  TASKTITLE='Generate a Free R set'
  DESCRIPTION="Generate a Free R set for a complete set of reflection indices to a given resolution (FreeRflag)"
  RANK = 2
  EXPORTMTZPARAMS = [ 'F_SIGF','FREEROUT' ]
  
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def drawContents(self):
    self.setProgramHelpFile('freerflag')
    self.openFolder(folderFunction='inputData')
    self.createLine(['advice', 'By default, this task will create a new set of freeR flags'])
    self.createLine(['widget', 'F_SIGF'])
    self.createLine(['widget', '-guiMode', 'radio', 'GEN_MODE'])
    self.createLine(['tip', 'FreeR column to be completed', 'widget','FREERFLAG'], toggle=['GEN_MODE', 'open', ['COMPLETE']])
    self.createLine(['label', 'Set fraction of reflections in freeR set', 'widget', 'FRAC',
                    'advice', 'Default fraction is 0.05',
                    'advice', 'Potential twinning operations will be taken into account'],
                    toggle=['GEN_MODE', 'open', ['GEN_NEW']])
    # self.openFolder(title='Advanced Options')

    # Maximum resolution options
    try:  # to cope with old jobs missing new parameters
      highResD = self.container.controlParameters.FSIGF_RESMAX.get()
      highResF = self.container.controlParameters.FREER_RESMAX.get()
      line = self.createLine(['advice',
                              'Maximum resolution (\xc5) in: data file ',
                              'label', 'data resolution',
                              'advice', ' FreeR file ',
                              'label', 'data resolution'],
                             toggleFunction=[self.getMaxReso,
                                             ['RELATIVE_RESOLUTION','GEN_MODE']])
      # get label widget, note itemAt counts from 0
      self.maximum_datares_label = line.layout().itemAt(1).widget()
      self.maximum_freeres_label = line.layout().itemAt(3).widget()
      self.createLine(['widget', 'CUTRESOLUTION',
                       'label','Cut resolution of FreeR to maximum resolution of data'],
                      toggle=['GEN_MODE', 'open', ['COMPLETE']])
    except:
      pass
    
    self.createLine(['label', 'Optionally, set high resolution limit at', 'widget', 'RESMAX', 'label', 'A.'], toggle=['UNIQUEIFY', 'open', [True]] )

    if self.isEditable():
      self.handle_gen_mode_changed()
      self.container.controlParameters.GEN_MODE.dataChanged.connect(self.handle_gen_mode_changed)
      qcb = getattr(self.getWidget('FREERFLAG'), 'jobCombo', None)
      if qcb: qcb.setItemText(qcb.findData(-1), '.. must be selected')
      self.container.inputData.FREERFLAG.dataChanged.connect( self.handle_freer_changed)

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  @QtCore.Slot()
  def handle_gen_mode_changed(self):
    #print(">> handle_gen_mode_changed", self.container.controlParameters.GEN_MODE)
    mode_new = self.container.controlParameters.GEN_MODE == 'GEN_NEW'
    self.container.inputData.FREERFLAG.setQualifier('allowUndefined', mode_new)
    self.getWidget('FREERFLAG').validate()
    self.handle_freer_changed()
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  @QtCore.Slot()
  def handle_freer_changed(self):
    highResD = -1.0
    highResF = -1.0
    
    if self.container.inputData.F_SIGF.fileContent.resolutionRange.high.isSet():
      highResD = float(self.container.inputData.F_SIGF.fileContent.resolutionRange.high)
    if self.container.inputData.FREERFLAG.fileContent.resolutionRange.high.isSet():
      highResF = float(self.container.inputData.FREERFLAG.fileContent.resolutionRange.high)

    # print("MaxReso:", highResD, highResF)
    self.container.controlParameters.FSIGF_RESMAX.set(highResD)
    self.container.controlParameters.FREER_RESMAX.set(highResF)

    self.container.controlParameters.CUTRESOLUTION.set(False)
    relativeResolution = 'UNKNOWN'  # Unknown
    if highResD > 0.0 and highResF > 0.0:
      if (highResD - highResF) > 0.002:
        print("FreeR goes to higher resolution than data")
        relativeResolution = 'FREERHIGHER'
        self.container.controlParameters.CUTRESOLUTION.set(True)
      elif (highResD - highResF) < 0.002:
        relativeResolution = 'DATAHIGHER'  # Data higher resolution than FreeR
      else:
        relativeResolution = 'SAME'
    self.container.controlParameters.RELATIVE_RESOLUTION.set(relativeResolution)

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  @QtCore.Slot()
  def getMaxReso(self):
    # return true if line should be visible
    if self.container.controlParameters.GEN_MODE == 'GEN_NEW':
      self.container.controlParameters.RELATIVE_RESOLUTION.set('UNKNOWN')
    if self.container.controlParameters.RELATIVE_RESOLUTION == 'UNKNOWN':
      if self.maximum_datares_label:
        self.maximum_datares_label.setText('')
        self.maximum_freeres_label.setText('')
      return False
    # Resolution values available
    highResD = self.container.controlParameters.FSIGF_RESMAX.get()
    highResF = self.container.controlParameters.FREER_RESMAX.get()
    self.maximum_datares_label.setText('{:5.2f}'.format(highResD))
    self.maximum_freeres_label.setText('{:5.2f}'.format(highResF))
    return True

