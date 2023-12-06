"""
     tasks/import_xia2.py: CCP4 GUI Project
     Copyright (C) 2013 STFC

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
     Liz Potterton Feb 2013 - Create basic import xia2 gui
"""

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets
 

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
class CTaskimport_xia2(CCP4TaskWidget.CTaskWidget):
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'import_xia2'
  TASKVERSION = 0.0
  #MN For now, I'd like to stick with the "AlternativeImportXIA2"
  #Thanks Liz for your patience with this decision
  #TASKMODULE='test'
  TASKTITLE='Import from Xia2 folder'
  DESCRIPTION = '''Import merged and unmerged X-ray reflections and results from Xia2'''
  CLONEABLE = False
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)


  def drawContents(self):

    self.setProgramHelpFile('import_xia2')
                        
    self.openFolder(folderFunction='inputData',title='Input Data', followFrom=False)

    self.createLine ( [  'advice' ,"Select either the 'top'  XIA2 directory to import all processing runs or select one processing run" ] )
    self.createLine ( [  'advice' ,'Each processing run in the XIA2 directory will be automatically imported as a separate job' ] )

    self.createLine ( [  'label' , 'XIA2 Directory',
                         'widget','XIA2_DIRECTORY' ] )

   
    
    #self.createLine ( [ 'widget',  'COPY_XIA2_DIRECTORY', 'label',
    #                    'Copy and save entire XIA2 directory in project?' ] )
