from ....qtgui.CCP4TaskWidget import CTaskWidget


class sheetbend_gui(CTaskWidget):
    TASKNAME = 'sheetbend' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'refinement' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='Shift field refinement'
    TASKTITLE='Fast preliminary refinement of atomic model coordinates or temperature factors, including at low resolution'
    DESCRIPTION = '''Fast preliminary refinement of atomic model using the program sheetbend'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['prosmart_refmac','coot_rebuild', 'modelcraft' ]

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
        
        self.createLine(['subtitle','Input reflections'])
        self.openSubFrame(frame=True)
        self.createLine( [ 'tip','Input reflections','widget','F_SIGF' ] )
        self.createLine( [ 'widget', 'FREERFLAG' ] )
        self.createLine( [ 'widget','XYZIN' ] )
        self.closeSubFrame()

        self.closeFolder ( )

        self.openFolder ( folderFunction='controlParameters', title='Basic Options' )

        self.createLine(['subtitle','Parameters to refine'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'widget', 'REFINE_COORD',   'label', 'coordinates' ] )
        self.createLine ( [ 'widget', 'REFINE_U_ISO',   'label', 'isotropic B factors' ] )
        self.createLine ( [ 'widget', 'REFINE_U_ANISO', 'label', 'anisotropic B factors' ] )
        self.closeSubFrame()
        self.createLine ( [ 'label', 'Number of cycles', 'widget', 'CYCLES' ] )
        self.createLine ( [ 'tip', 'Enter a single number or a comma separated list', 'label', 'Resolution', 'widget', 'RESOLUTION', 'label', 'Angstroms' ] )

        self.closeFolder ( )

        self.openFolder ( folderFunction='controlParameters', title='Advanced Options' )

        self.createLine ( [ 'label', 'Sphere radius', 'widget', 'RADIUS_SCALE', 'label', 'times resolution' ] )

        self.createLine(['subtitle','Refine-regularize macro cycles'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'label', 'Perform pseudo-regularization', 'widget', 'PSEUDO_REGULARIZE' ] )
        self.createLine ( [ 'label', 'Number of refine-regularize macro-cycles', 'widget', 'REFINE_REGULARIZE_CYCLES' ] )
        self.createLine(['subtitle','Additional parameters to refine before regularization'])
        self.createLine ( [ 'widget', 'POSTREFINE_COORD',   'label', 'coordinates' ] )
        self.createLine ( [ 'widget', 'POSTREFINE_U_ISO',   'label', 'isotropic B factors' ] )
        self.createLine ( [ 'widget', 'POSTREFINE_U_ANISO', 'label', 'anisotropic B factors' ] )
        self.closeSubFrame()

        self.closeFolder ( )
