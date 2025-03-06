"""
     tasks/privateer/privateer_gui.py: CCP4 GUI Project
     Copyright (C) 2014 University of York

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
     Jon Agirre         2014 - Started development
     Jon Agirre         2016 - Revamped for 7.0 release
     Jon Agirre         2019 - New version to support MKIV functionality
     Haroldas Bagdonas  2021 - New version Privateer MKIV with new functionality.

"""
from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class privateer_gui(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'privateer'
  TASKVERSION = 0.1
  TASKMODULE='validation'
  TASKTITLE='Validation of carbohydrate structures - Privateer'
  DESCRIPTION='Validation, re-refinement and graphical analysis of carbohydrate structures'
  SHORTTASKTITLE='Privateer'

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('privateer')


# the input data tab starts here

    folder = self.openFolder(folderFunction='inputData',title='Input Data')

    self.createLine( [ 'subtitle', 'Model and experimental data', 'Your model must contain monosaccharides. Please check the documentation for a list of known monosaccharide three-letter codes. Reflection data are required for calculating correlation between model and unbiased density.' ] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'widget', 'XYZIN'])
    self.createLine( [ 'widget', 'F_SIGF' ] )
    self.createLine( [ 'label', 'Change the mask radius around the sugar atoms to', 'widget', 'RADIUSIN', 'label', 'Angstroems' ] )
    # self.createLine( [ 'widget', 'BLOBS', 'label', 'EXPERIMENTAL FEATURE: Scan for unmodelled Glycosylation, change Electron Density Blob Threshold level to', 'widget', 'BLOBSLEVEL', 'label' ])
    
    self.createLine( [ 'label', '&nbsp;' ])
    self.createLine( [ 'widget', 'NEW_SUGAR', 'label', 'A sugar I want to validate is not yet part of the Chemical Component Dictionary'])
    self.createLine( [ 'subtitle', 'Define a new sugar', 'Please define your sugar using the fields below. If your sugar is a polysaccharide, you will have to enter information corresponding to just one ring.' ], toggle=[ 'NEW_SUGAR', 'open', [True]] )
    self.openSubFrame ( frame=[True], toggle=[ 'NEW_SUGAR', 'open', [True]] )
    self.createLine( [ 'label', 'Analyse sugar with code', 'widget', 'CODEIN', 'label', 'and type', 'widget', 'ANOMER', 'widget', 'HAND','widget', 'RING_TYPE' ] )
    self.createLine( [ 'label', 'Expected minimal energy ring conformation:', 'widget', 'CONFORMATION_PYRANOSE' ], toggle=[ 'RING_TYPE', 'open', 'pyranose'] )
    self.createLine( [ 'label', 'Expected minimal energy ring conformation:', 'widget', 'CONFORMATION_FURANOSE' ], toggle=[ 'RING_TYPE', 'open', 'furanose'] )
    self.createLine( [ 'label', 'Ring atoms: ', 'widget', 'RING_OXYGEN', 'widget', 'RING_C1', 'widget', 'RING_C2', 'widget', 'RING_C3', 'widget', 'RING_C4', 'widget', 'RING_C5' ], toggle=[ 'RING_TYPE', 'open', 'pyranose'] )
    self.createLine( [ 'label', 'Ring atoms: ', 'widget', 'RING_OXYGEN', 'widget', 'RING_C1', 'widget', 'RING_C2', 'widget', 'RING_C3', 'widget', 'RING_C4' ], toggle=[ 'RING_TYPE', 'open', 'furanose'] )
    self.closeSubFrame ( )

    self.container.controlParameters.NEW_SUGAR.dataChanged.connect(self.updateRequirements )
    self.closeSubFrame ( )
    self.closeFolder()

    folder = self.openFolder(title='Settings')
    self.createLine( [ 'subtitle', 'Glycosylation analysis', 'Here you can set up how the diagrams will look like. The schemes will follow the Essentials of glycobiology 3rd edition notation with a choice of colours. They are vector-based and can be used in publications straight away.' ] )
    self.openSubFrame ( frame=[True] )

    self.createLine( [ 'label', 'Create glycan diagrams in', 'widget', 'OLDSTYLE', 'label', 'style, in', 'widget', 'VERTICAL', 'label', 'orientation' ])
    self.createLine( [ 'label', 'Colour them using the' , 'widget', 'ESSENTIALS' , 'label' , 'colour scheme with ', 'widget', 'INVERT', 'label', 'outlines' ])
    self.createLine( [ 'label', 'Validate glycan structures assuming a', 'widget', 'EXPRESSION', 'label', 'expression system' ])
    self.createLine( [ 'widget', 'GLYTOUCAN', 'label', 'Validate Glycans with GlyTouCan and GlyConnect databases'])
    self.openSubFrame ( frame=[False], toggle=[ 'GLYTOUCAN', 'open', [True]] )
    self.createLine( [ 'widget', 'CLOSESTMATCH', 'label', 'GlyConnect: Don\'t look for the closest match, if modelled glycan is not found in the database.'])
    self.createLine( [ 'widget', 'ALLPERMUTATIONS', 'label', 'Conduct all possible permutation combinations(WARNING: Will take a long time to finish, should only be really used for O-Glycans).'])
    self.closeSubFrame ( )

    self.createLine( [ 'label', '&nbsp;' ])
    self.createLine( [ 'subtitle', 'Performance settings', 'Here you can tweak performance related settings, such as whether to run Privateer with a single CPU thread or customize the number of CPU threads to run Privateer with. Privateer is tweaked by default to try to obtain best performance possible that is offered by the computer, i.e. use maximum number of threads available.' ] )
    self.openSubFrame ( frame=[True] )
    self.createLine( [ 'label', 'Run Privateer with', 'widget', 'NUMTHREADS', 'label', 'threads.' ] )
    self.closeSubFrame ( )

    self.closeFolder()
    # self.closeFolder()

  @QtCore.Slot()
  def updateRequirements ( self ) :
    if self.container.controlParameters.NEW_SUGAR :
        self.container.controlParameters.CODEIN.setQualifier ( 'allowUndefined', False )
        self.container.controlParameters.CODEIN.validate()
    else :
        self.container.controlParameters.CODEIN.setQualifier ( 'allowUndefined', True )
        self.container.controlParameters.CODEIN.validate()
