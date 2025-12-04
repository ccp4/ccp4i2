"""
    ctruncate_gui.py
    Copyright (C) 2015STFC
    Author: Liz POtterton

"""

from baselayer import QtGui, QtWidgets,QtCore
from qtgui.CCP4TaskWidget import CTaskWidget
from core.CCP4ErrorHandling import CErrorReport

class ctruncate_gui(CTaskWidget):

    TASKNAME = 'ctruncate'
    TASKVERSION = 0.1
    TASKMODULE='expt_data_utility'
    TASKTITLE='Convert intensities to amplitudes'
    DESCRIPTION = 'Convert reflection intensities to structure factors (ctruncate)'

    ERROR_CODES = { 200 : { 'description' : 'Selected input data is not intensities' }
                    }

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):

        self.setProgramHelpFile ( 'ctruncate' )
        
        self.openFolder ( folderFunction='inputData',followFrom=False )

        self.container.inputData.OBSIN.setQualifier('allowUndefined',False)
        self.container.inputData.HKLIN.setQualifier('allowUndefined',True)
        self.createLine ( [ 'subtitle', 'Reflection data - mean intensities or anomalous intensities','Select reflection data or import an MTZ file'] )
        self.createLine ( [ 'widget', 'OBSIN' ] )

        self.openFolder ( title='Options' )
        
        self.createLine ( [ 'subtitle', 'AU content as sequence or number of residues','Scaling is improved by knowing AU content'] )
        self.openSubFrame(frame=[True])
        self.createLine ( [ 'widget', 'SEQIN' ] )
        self.createLine ( [ 'label' , 'OR - number of residues in the asymmetric unit', 'widget', 'NRES' ] )
        self.closeSubFrame()
        
        self.createLine ( ['widget', 'NO_ANISO', 'label', 'Do not perform anisotropy correction' ] )
# AL: the next line is not functional, Imean is available from aimless, therefore commented out:
#       self.createLine ( ['widget', 'OUTPUT_INTENSITIES', 'label', 'Output Imean for anomalous data'] )

    def fix(self):
      self.container.inputData.OBSIN.setContentFlag()
      if self.container.inputData.OBSIN.contentFlag.__int__() not in [ self.container.inputData.OBSIN.CONTENT_FLAG_IPAIR, self.container.inputData.OBSIN.CONTENT_FLAG_IMEAN]:
          return CErrorReport(self.__class__,200,stack=False)
      self.container.controlParameters.OUTPUTMINIMTZ.set(True)
      self.container.inputData.OBSIN.setContentFlag()
      if  self.container.inputData.OBSIN.contentFlag==1:
        self.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG.set(2)
      else:
        self.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG.set(4)
      return CErrorReport()
