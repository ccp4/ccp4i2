"""
    AcedrgLink_gui.py: CCP4 GUI Project
    
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

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class AcedrgLink_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'AcedrgLink' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    SHORTTASKTITLE='AceDRG in link generation mode'
    TASKTITLE='AceDRG in link generation mode'
    DESCRIPTION = '''AceDRG in link generation mode'''
    WHATNEXT = ['coot_rebuild']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
