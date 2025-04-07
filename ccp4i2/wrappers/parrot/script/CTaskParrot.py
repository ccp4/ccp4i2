"""
Copyright (C) 2012 STFC
Liz Potterton Aug 2012 = Parrot gui
"""

from ....core import CCP4ErrorHandling
from ....qtgui import CCP4TaskWidget


class CTaskParrot(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'parrot'
  TASKVERSION = 0.0
  TASKMODULE= 'density_modification' 
  TASKTITLE='Density modification - PARROT'
  SHORTTASKTITLE='PARROT'
  DESCRIPTION='Modify the electron density (Parrot)'
  WHATNEXT = [ 'coot_rebuild',['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml']]
  MGDISPLAYFILES = ['FPHIOUT']
  EXPORTMTZPARAMS = ['F_SIGF','ABCD','ABCDOUT']
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('parrot')

    folder = self.openFolder(folderFunction='inputData',title='Input Data')
    
    self.createLine(  [ 'subtitle', 'Select experimental data', 'Observed structure factors and initial phasing (e.g. from experimental phasing or molecular replacement) are required' ] )

    self.openSubFrame( frame=[True] )
    self.createLine(  [ 'widget','F_SIGF' ] )
    self.createLine(  [ 'widget','ABCD' ] )
    self.closeSubFrame()

    self.createLine(  [ 'subtitle', 'Select AU content for solvent content estimation', 'If you have additional information about the solvent content, you can override this under basic options'] )
    self.openSubFrame( frame=[True] )
    self.createLine(  [ 'widget','ASUIN' ] )
    self.closeSubFrame()

    self.createLine(  [ 'subtitle', 'Select NCS information', 'Provide a heavy atom or partial model to allow NCS determination and averaging' ] )
    self.openSubFrame( frame=[True] )
    self.createLine( ['tip','Provide a heavy atom or partial model to allow NCS determination and averaging','widget','-guiMode','radio','XYZIN_MODE',] )
    stack = self.openStack(controlVar='XYZIN_MODE')
    self.createLine(['label','No NCS model'])
    self.createLine(['widget','XYZIN_HA'])
    self.createLine(['widget','XYZIN_MR'])
    self.closeStack()
    self.closeSubFrame()


    folder = self.openFolder(folderFunction='controlParameters',title='Basic Options')
    self.createLine(['label', 'Number of cycles', 'widget', 'CYCLES'])
    self.createLine(['label','Override fractional solvent content','widget','SOLVENT_CONTENT','advice','  Determined from sequence if blank'])

    folder = self.openFolder(folderFunction='controlParameters',title='Advanced Options')

    self.createLine(  [ 'subtitle', 'Free-R flag' ] )
    self.createLine(['widget','FREERFLAG'])
    self.createLine(  [ 'subtitle', 'Map coefficients for starting map' ] )
    self.createLine(['widget','F_PHI'])

    self.createLine(  [ 'subtitle', 'Additional options' ] )
    self.createLine(['widget','ANISOTROPY_CORRECTION','label','Apply anisotropy correction'])

    self.createLine(['label','Truncate data beyond resolution limit','widget','RESOLUTION','label','(Angstroms)','advice','  Use all data if blank'])
    self.createLine(['label','Radius of density sphere for testing NCS operators','widget','NCS_MASK_FILTER_RADIUS','label','(Angstroms)'])
        
    folder = self.openFolder(folderFunction='controlParameters',title='Reference structures ')
    self.createLine( [ 'advice', 'You should normally let Parrot choose reference structures'] )
    self.createLine( [ 'widget','F_SIGF_REF' ] )
    self.createLine( [ 'widget','ABCD_REF' ] )
    self.createLine( [ 'widget','XYZIN_REF' ] )

    #for w in self.findChildren(CCP4TaskWidget.CStackedWidget):
    #  print '   ',w,w.controlVar

    #self.updateViewFromModel()
    
  def fix(self):
    # Unset possible input coord files that will not be used unless in appropriate
    # NCS mode
    inp = self.container.inputData

    if inp.XYZIN_HA.isSet() and inp.XYZIN_MODE != 'ha': inp.XYZIN_HA.unSet()
    if inp.XYZIN_MR.isSet() and inp.XYZIN_MODE != 'mr': inp.XYZIN_MR.unSet()

    return CCP4ErrorHandling.CErrorReport()
