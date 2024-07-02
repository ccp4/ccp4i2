"""
     tasks/pyphaser_mr/CTaskPyphaser_mr.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from qtgui.CCP4TaskWidget import CTaskWidget
from PySide6 import QtCore

#-------------------------------------------------------------------
class CTaskPyphaser_mr(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'pyphaser_mr'
  TASKVERSION = 0.1
  TASKMODULE='test'
  TASKTITLE='MR using Phaser (pythonic)'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def setDefaultParameters(self):
    # Reimplement in tasks to set initial parameters

    # read ccp4/src/phaser/source/phaser/defaults for defaults

    pass

  def drawContents(self):

    self.setProgramHelpFile('phaser')

    self.openFolder(folderFunction='inputData')
     
    self.openSubFrame()
    self.createLine ( [ 'widget', 'F_SIGF' ] )
    self.closeSubFrame()

    
    self.setMenuText('COMP_BY',{ 'DEFAULT':'Provide no guidance to Phaser',
                     'MW':'Provide an estimate of the molecular weight of protein and nucleic acid in ASU',
                     'ASU':'Provide a full specification of the ASU content'
                     })
    self.createLine ( ['label', 'For estimating asymmetric unit contents:','widget','COMP_BY'])
    self.createLine ( ['label', '\n '], toggle=['COMP_BY','open',['DEFAULT']])
                     
    self.createLine ( [ 'widget', '-title','Contents of asymmetric unit', 'ASU_COMPONENTS' ], toggle=['COMP_BY','open',['ASU']])
    self.createLine ( [ 'label','Molecular weight (Da) of protein in the ASU','stretch','widget', 'ASU_PROTEIN_MW' ], toggle=['COMP_BY','open',['MW']])
    self.createLine ( [ 'label','Molecular weight (Da) of nucleic acid','stretch','widget', 'ASU_NUCLEICACID_MW' ], toggle=['COMP_BY','open',['MW']])
    #The validity of the ASU_COMPONENTS widget contents are dependent on COMP_BY
    self.container.controlParameters.COMP_BY.dataChanged.connect(self.getWidget('ASU_COMPONENTS').validate)
        
    self.createLine ( [ 'label', 'Use the ensembles below as', 'widget' ,  '-guiMode', 'radio', 'SEARCHMODE' ] )
    self.createLine ( [  'widget', '-title','List of search model ensembles','ENSEMBLES' ] )

    
    self.openFolder(title='Important Options')

    self.setMenuText('MODE',{ 'MR_AUTO':'fully automated',
                              'MR_FRF':'rotation function only',
                              'MR_FTF':'translation function only'})
    self.createLine( [ 'label', 'MR method',
                       'widget', 'MODE',
                       ] )

    self.setMenuText('SGALT_SELECT',{ 'NONE':'Use input spacegroup only',
                 'HAND':'Test spacegroup and enantiomorph',
                 'ALL':'Test all spacegroups in point group',
                 'LIST':'Select from list'})
    self.createLine( [ 'label', 'Choice of spacegroups to test',
                                   'widget', 'SGALT_SELECT',
                                   ] )
    self.setToolTip('SGALT_TEST','Enter a list of space-separated spacegroups')
    self.createLine( [ 'widget', 'SGALT_TEST' ], toggle=['SGALT_SELECT','open',['LIST']])
    #The validity of the SGALT_TEST widget contents are dependent on SGALT_SELECT
    self.container.controlParameters.SGALT_SELECT.dataChanged.connect(self.getWidget('SGALT_TEST').validate)
    
    self.openFolder(title='Additional Options')

    self.setToolTip('RESOLUTION_HIGH','Limiting the resolution will make the job faster but less sensitive')
    self.createLine ( [ 'label' , 'Set high resolution limit to ',
                        'widget', 'RESOLUTION_HIGH' ] )

    self.setToolTip('PEAKS_ROT_CUTOFF','Lower this to try more rotation solutions (slower)')
    self.createLine ( [ 'label' , 'Select rotation peaks over',
                        'widget', 'PEAKS_ROT_CUTOFF',
                        'label', 'percent.' ] )

    self.setToolTip('PACK_CUTOFF','Increase this if inaccurate loops are causing packing test to fail')
    self.createLine ( [ 'label' , 'Maximum clashes allowed as percentage of total Calpha',
                        'widget', 'PACK_CUTOFF' ] )

    self.openFolder(title='Advanced Options')

    self.createLine ( [ 'advice', 'Use more local processors if available' ] )
    self.createLine ( [ 'label', 'Number of processors',
                        'widget', 'NJOBS' ] )

    self.createLine ( [ 'label' , 'Number of top solutions to output',
                        'widget', 'NUM_SOL_OUT' ] )

    self.createLine ( [ 'widget', 'PERMUTATIONS',
                        'label', 'Permute search set' ] )
  
   
