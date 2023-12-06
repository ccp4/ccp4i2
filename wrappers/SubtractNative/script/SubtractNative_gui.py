"""
    SubtractNative_gui.py: CCP4 GUI Project
    
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
class SubtractNative_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'SubtractNative' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'refinement' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='Subtract native map'
    TASKTITLE='Subtract a defied fraction of an atom map from an input map'
    DESCRIPTION = '''Reads a set of map coefficients, generates a corresponding map, and modifies it by subtracting some fraction of the FC map corresponding to a set of coordinates that have been provided. The output is a set of modified map coefficients'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
        
        self.openSubFrame(frame=True, title='Inputs')
        self.createLine ( [ 'tip','Input map coefficients','widget','MAPIN' ] )
        self.createLine ( [ 'tip','Input model','widget','XYZIN' ] )
        self.getWidget('XYZIN').showAtomSelection()
        self.createLine ( [ 'tip','Fractionof FC to subtract','widget','FRACTION' ] )
        self.closeSubFrame()
        
