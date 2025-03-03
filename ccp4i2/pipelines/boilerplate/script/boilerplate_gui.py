"""
    ZZPipelineNameZZ_gui.py: CCP4 GUI Project
    
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
class ZZPipelineNameZZ_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'ZZPipelineNameZZ' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'ZZFolderZZ' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='ZZShortTitleZZ'
    TASKTITLE='ZZLongTitleZZ'
    DESCRIPTION = '''ZZDescriptionZZ'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
        
        self.createLine(['subtitle','Input coordinates'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input model','widget','XYZIN_LIST' ] )
        #self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()
        
        self.createLine(['subtitle','Input reflections'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input reflections','widget','F_SIGF' ] )
        self.closeSubFrame()
