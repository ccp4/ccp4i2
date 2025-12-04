from __future__ import print_function


"""
     tasks/molrep_mr/Cmolrep_mr.py: CCP4 GUI Project
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
     Andrey Lebedev September 2011 - molrep_mr gui
"""

from baselayer import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets



class Cmolrep_mr(CCP4TaskWidget.CTaskWidget):

  TASKNAME = 'molrep_mr'
  TASKVERSION = 0.0
  TASKMODULE='test'
  WHATNEXT = ['prosmart_refmac','coot_rebuild']
  TASKTITLE = 'Molecular Replacement and refinement- MOLREP'
  SHORTTASKTITLE='MOLREP MR'
  DESCRIPTION='Molecular replacement (Molrep)'
  MGDISPLAYFILES = ['XYZOUT']  
  PROGRAMHELP = ['molrep']
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)


  def drawContents(self):

    self.setProgramHelpFile('molrep_mr')

#-  --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='inputData',title='Input Data')


    """
    self.createLine( [ 'advice', 'What to do?' ] )
    self.setMenuText( 'PERFORM', {
       'srf': 'Self Rotation Function',
       'pat': 'Rotation and Translation Searches',
       'den': 'Search in Density after Refinement',
    } )
    self.createLine( [ 'widget', '-guiMode', 'multiLineRadio', 'PERFORM' ] )
    """
    

    self.openSubFrame( toggle = [ 'PERFORM','close', [ 'den' ]] )
    self.createLine( [ 'advice', 'Experimental Data' ] )
    self.createLine( [ 'widget', 'F_SIGF' ] )
    self.closeSubFrame()

    self.openSubFrame( toggle=[ 'PERFORM','close', [ 'srf' ]] )
    self.createLine( [ 'advice', 'Search Model' ] )
    self.createLine( [ 'widget', 'XYZIN' ] )
    self.createLine( [ 'widget', 'SEQIN' ] )
    self.createLine( [ 'label', 'The number of monomers to search for', 'widget', 'NMON' ] )
    print('molrep_mr setWhatsThis',self.getWidget('NMON'))
    self.getWidget('NMON').setWhatsThis("The number of monomers in the asymmetric unit - recommended that you choose 'auto' and let the program decide.")
    self.closeSubFrame()

    self.openSubFrame( toggle = [ 'PERFORM','open', [ 'pat' ]] )
    self.createLine( [ 'advice', 'Fixed Model' ] )
    self.createLine( [ 'widget', 'XYZIN_FIX' ] )
    self.closeSubFrame()

    '''
    self.openSubFrame( toggle = [ 'PERFORM','open', [ 'den' ] ])
    self.createLine( [ 'advice', 'Partial Model and Map Coefficients: define PDB and MTZ after the same refmac run' ] )
    self.createLine( [ 'widget', 'XYZIN_FIX' ] )
    self.createLine( [ 'widget', 'F_PHI_MAP' ] )
    self.closeSubFrame()
    '''

#-   --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='controlParameters',title='Basic Options')

    self.createLine ( [ 'advice','Searching Laue group currently does not output a model - consult log file' ])

    self.createLine ( [ 'label','Search in space group(s)','widget', 'SG_OPTIONS'] )
    
    self.createLine ( [ 'widget','-label','Search in alternative space group','SG'], toggle=['SG_OPTIONS','open',[ 'specify' ] ]  )
    
    self.createLine( [ 'label', 'Use data in resolution range from low ', 'widget', 'RESMIN', 'label', 'to high','widget', 'RESMAX' ] )

    self.createLine( [ 'advice', 'Perform alignment and use it to rename residues and trim side chains' ] )
    self.setMenuText( 'SEQ', {
       'y': 'always',
       'd': 'only for sequence identity > 20%',
       'n': 'never'
    } )
    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'SEQ' ] )


    self.createLine( [ 'advice', 'B-factors modification options' ] )
    self.setMenuText( 'SURF', {
       'y': 'Increase B-factor on the molecular surface for all functions',
       'c': 'Increase B-factor on the surface for Packing function only',
       'n': 'Do not do anything',
       '2': 'Set B-factors of all atoms to 20',
       'a': 'Use poly-alanine model with all B-factors 20',
    } )
    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'SURF' ], toggle=[ 'OPEN_SURF', 'open', [ True ] ] )

    self.createLine( [ 'advice', 'Number of peaks to analyse' ] )
    self.createLine( [ 'label', '         ', 'label', 'Number of Rotation Function peaks', 'widget', 'NP'] )
    self.createLine( [ 'label', '         ', 'label', 'Number of Translation Function peaks', 'widget', 'NPT'])


    folder = self.openFolder(folderFunction='controlParameters',title='Advanced Options')

    self.openSubFrame( toggle = [ 'PERFORM','open', [ 'den' ]] )
    self.createLine( [ 'advice', 'Search method' ] )
    self.createLine( [ 'label', '         ', 'label', 'SAPTF = Spherically Averaged Phased Translation Function' ],)
    self.createLine( [ 'label', '         ', 'label', 'PRF = Phased Rotation Function' ] )
    self.createLine( [ 'label', '         ', 'label', 'RF(M) = Rotation Function (f-obs from the density outside of the fixed Model)' ] )
    self.createLine( [ 'label', '         ', 'label', 'RF(S) = Rotation Function (f-obs from the density inside a Sphere)' ] )
    self.setMenuText( 'PRF', {
       'n': 'RF(M) + PTF',
       'y': 'SAPTF + PRF + PTF',
       's': 'SAPTF + RF(S) + PTF',
    } )
    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'PRF' ] )
    self.closeSubFrame()


    self.openSubFrame( toggle = [ 'PERFORM','open', [ 'pat' ] ] )
    self.createLine( [ 'advice', 'Scoring putative solutions' ] )
    self.createLine( [ 'label', '         ', 'label', 'CC = Correlation Coefficient' ] )
    self.createLine( [ 'label', '         ', 'label', 'PF = Packing Function' ] )
    self.setMenuText( 'SCORE', {
       'y': 'Use CC times PF as a score and stop translation search if contrast is > 3.0',
       'n': 'Use CC times PF and do not stop',
       'c': 'Use CC and do not stop',
    } )
    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'SCORE' ] )
    self.createLine( [ 'label', '         ', 'label', 'Expected number of copies (for contrast calculation only)', 'widget', 'NMON_EXP'] )


    self.createLine( [  'advice', 'Scaling' ] )
    self.setMenuText( 'ANISO', {
       'y': 'anisotropic',
       'n': 'isotropic',
       'k': 'none',
    } )
    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'ANISO' ] )
    self.closeSubFrame()


    self.createLine( [ 'advice', 'High pass filter parameter (B-add, the B-factor applied to input structure amplitudes)' ] )

    self.setMenuText( 'HIGH_PATH_VAR', {
       's': 'From identity between model and sequence (if sequence given)',
       'i': 'From identity specified manually',
       'r': 'From high resolution limit',
       'b': 'Directly as the value of additional B-factor',
    } )
    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'HIGH_PATH_VAR' ])
    lab1 = 'Identity between search and target sequences (from 0 to 1)'
    self.createLine( [ 'label', '         ', 'label', lab1, 'widget', 'SIM' ], toggle=['HIGH_PATH_VAR', 'open' , ['i'] ]  )
    lab1 = 'High resolution limit, Ang'
    self.createLine( [ 'label', '         ', 'label', lab1, 'widget', 'RESMAX' ], toggle= [ 'HIGH_PATH_VAR', 'open' , ['r'] ] )
    lab1 = 'B-add'
    self.createLine( [ 'label', '         ', 'label', lab1, 'widget', 'BADD' ], toggle=[ 'HIGH_PATH_VAR', 'open' , ['b'] ] )


    self.createLine( [ 'advice', 'Low pass filter parameter (B-off, the B-factor of the removed fraction of structure amplitudes)' ] )

    self.setMenuText( 'LOW_PATH_VAR', {
       'c': 'From completeness of the search model',
       'r': 'From low resolution limit',
       'b': 'Directly as the value of additional B-factor',
    } )
    self.createLine( [ 'label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'LOW_PATH_VAR' ])



    '''
    print 'CTaskMolrep stackedWidgets'
    for w in self.findChildren( CCP4TaskWidget.CStackedWidget ) :
      print '   ', w, w.controlVar

    self.updateViewFromModel()
    '''

    self.container.guiParameters.OPEN_HIGH_PATH_VAR = True



