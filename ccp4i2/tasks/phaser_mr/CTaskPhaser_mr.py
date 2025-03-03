from __future__ import print_function

"""
     tasks/phaser_mr/CTaskPhaser_mr.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from qtgui.CCP4TaskWidget import CTaskWidget
from PySide2 import QtCore

#-------------------------------------------------------------------
class CTaskPhaser_mr(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'phaser_mr'
  TASKVERSION = 0.1
  TASKMODULE='test'
  TASKTITLE='MR using Phaser'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def drawContents(self):
    print('Whoops')

    self.setProgramHelpFile('phaser')

    self.openFolder(folderFunction='inputData')
     
    self.openSubFrame()
    self.createLine ( ['label', '\n '])
    self.createLine ( [ 'widget', 'F_SIGF' ] )
    self.closeSubFrame()

    self.setMenuText('COMP_BY',{ 'DEFAULT':'Provide no guidance to Phaser',
                 'MW':'Provide an estimate of the molecular weight of protein and nucleic acid in ASU',
                 'ASU':'Provide a full specification of the ASU content by sequence'
                 })
    self.createLine(['advice',''])
    self.createLine(['advice',''])
    self.createLine(['advice',"Phaser's maximum likelihood estimates will be best if you provide guidance about the contents of the asymmetric unit"])
    self.createLine(['advice',"Phaser also uses this information  to carry out its AU content analysis, to suggest how many copies to expect"])
    self.createLine ( ['label', 'For estimating asymmetric unit contents:','widget','COMP_BY'])
    self.createLine ( ['label', '\n '], toggle=['COMP_BY','open',['DEFAULT']])
    
    self.createLine ( [ 'widget', '-title','Contents of asymmetric unit  - click "Show list" if more than one type of chain is present', 'ASU_COMPONENTS' ], toggle=['COMP_BY','open',['ASU']])
    self.createLine ( [ 'label','Molecular weight (Da) of protein in the ASU','stretch','widget', 'ASU_PROTEIN_MW' ], toggle=['COMP_BY','open',['MW']])
    self.createLine ( [ 'label','Molecular weight (Da) of nucleic acid','stretch','widget', 'ASU_NUCLEICACID_MW' ], toggle=['COMP_BY','open',['MW']])
    #The validity of the ASU_COMPONENTS widget contents are dependent on COMP_BY
    self.container.controlParameters.COMP_BY.dataChanged.connect(self.getWidget('ASU_COMPONENTS').validate)
    
    self.createLine(['advice',''])


    self.createLine ( [  'widget', '-title','Search model(s) - click "Show list" if more than one copy or more than one search model','ENSEMBLES' ] )
    self.createLine ( [ 'label', 'Use the ensembles above as', 'widget' ,  '-guiMode', 'radio', 'SEARCHMODE' ] )
    
    self.createLine(['advice',''])
    self.createLine(['advice',''])
    self.createLine(['advice','If some part of the structure has already been positioned in the unit cell, you can enter it in the "Known structure" tab'])

    self.openFolder(title='Known structure')
    self.createLine ( ['label', '\n '])
    self.createLine ( [  'widget', '-title','Already positioned structure - click "Show list" if more than one copy or more than one search model','FIXED_STRUCTURE' ] )

    
    self.openFolder(title='Main Options')

    self.setMenuText('MODE',{ 'MR_AUTO':'fully automated',
                              'MR_FRF':'rotation function only',
                              'MR_FTF':'translation function only',
                              'MR_RNP':'refine ensembles in a solution and phase accounting for anisotropy and tNCS',
                              'MR_RGR':'refine chains in an ensemble of a solution and phase accounting for anisotropy and tNCS'
                     })
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
    self.createLine( [ 'widget', 'SGALT_TEST' ],
        toggle=['SGALT_SELECT','open',['LIST']])
    #The validity of the SGALT_TEST widget contents are dependent on SGALT_SELECT
    self.container.controlParameters.SGALT_SELECT.dataChanged.connect(self.getWidget('SGALT_TEST').validate)

    
    self.openFolder(title='Extra Options')

    self.setToolTip('RESOLUTION_HIGH','Limiting the resolution will make the job faster but less sensitive')
    self.createLine ( [ 'label' , 'Set high resolution limit to ',
                        'widget', 'RESOLUTION_HIGH' ] )
    self.createLine ( [ 'label' , 'Set high resolution limit for final stage of AUTO_MR',
                       'widget', 'RESOLUTION_AUTO_HIGH' ] )

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
  
  
  def isValid(self):
      #Here override logic of whether this is a valid task to allow for CSeqDataFile from the
      #CASUComponentList being required ONLY IF COMP_BY has the value "ASU"
      invalidElements = CTaskWidget.isValid(self)
      from core import CCP4ModelData, CCP4XtalData
      widgLib = {"COMP_BY":"Not set yet"}
      self.getParams(widgLib)
      if widgLib["COMP_BY"] != "ASU":
          invalidElements = [invalidElement for invalidElement in invalidElements if (type(invalidElement) != CCP4XtalData.CAsuComponent and type(invalidElement) != CCP4ModelData.CSeqDataFile)]
      return invalidElements



   
