from ....core import CCP4Container
from ....qtgui.CCP4TaskWidget import CTaskWidget


#-------------------------------------------------------------------
class phaser_phil_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'phaser_phil' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'test' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='phaser_phil'
    TASKTITLE='Phaser auto generated GUI'
    DESCRIPTION = '''Phaser auto generated GUI'''
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
        self.closeSubFrame()

        # Basic parameters
        self.createLine(['tip','Subset of Phaser control parameters','subtitle',
            'Basic parameters'])
        self.nestedAutoGenerate(container=self.container.controlParameters,
            expertLevel=['0'])

        # Advanced parameters in a separate folder
        folder = self.openFolder(folderFunction='controlParameters',
            title='Advanced parameters',drawFolder=self.drawAdvanced)

    # Copied from xia2_dials_gui - if useful, should be moved to the shared
    # base class
    def drawAdvanced(self):
      self.nestedAutoGenerate(container=self.container.controlParameters,
            expertLevel=['0', '1'])

    # Copied from xia2_dials_gui - if useful, should be moved to the shared
    # base class
    def nestedAutoGenerate(self, container, expertLevel=['0']):
        """Autogenerate a GUI from nested containers using only elements at
        the specified expertLevels"""

        dataOrder = container.dataOrder()
        contents = [getattr(container,name) for name in dataOrder]

        # Filter by expertLevel
        expert = [c.qualifiers('guiDefinition').get('expertLevel', '0') \
                  for c in contents]
        dataOrder = [d for d, e in zip(dataOrder, expert) if e in expertLevel]
        contents  = [c for c, e in zip(contents, expert) if e in expertLevel]

        # If this container has any content other than nested containers,
        # autogenerate that content
        data = [name for name, model in zip(dataOrder, contents) \
                if not isinstance(model, CCP4Container.CContainer)]
        subcontainers = [model for model in contents \
                         if isinstance(model, CCP4Container.CContainer)]
        if data:
            # only open a sub-frame if there is a label available
            make_subframe = container.qualifiers('guiLabel') is not NotImplemented
            if make_subframe: self.openSubFrame(frame=True)
            self.autoGenerate(container, selection={'includeParameters':data})
            for c in subcontainers:
                self.nestedAutoGenerate(c, expertLevel=expertLevel)
            if make_subframe: self.closeSubFrame()
        else:
            # just drill down without opening a sub frame
            for c in subcontainers: self.nestedAutoGenerate(c,
                expertLevel=expertLevel)
        return
