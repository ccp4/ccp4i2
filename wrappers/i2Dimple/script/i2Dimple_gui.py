"""
    i2Dimple_gui.py: CCP4 GUI Project
    
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
class i2Dimple_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'i2Dimple' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'molecular_replacement' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    TASKTITLE='DIMPLE - Simple reflections + coordinates to map pipeline'
    SHORTTASKTITLE='DIMPLE'
    DESCRIPTION = '''This task uses the DIMPLE pipeline to generate maps for a new dataset, provided a possible starting model is available.  For simple cases of isomorphous data, the pipeline will use rigid body refinement in REFMAC to 'tweak' the starting model. Where unit cells are incompatible, it will attempt automated molecular replacement.'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
        self.createLine(['subtitle','Input coordinates'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input model','widget','XYZIN' ] )
        self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()
        self.createLine(['subtitle','Input reflections'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input reflections','widget','F_SIGF' ] )
        self.createLine ( [ 'tip','Free R','widget','FREERFLAG' ] )
        self.closeSubFrame()
        self.createLine(['subtitle','Control parameters'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'label','R-factor threshold to force MR','widget','MR_WHEN_R' ] )
        self.closeSubFrame()
