from ....qtgui.CCP4TaskWidget import CTaskWidget


class zanuda_gui(CTaskWidget):
    TASKMODULE = 'refinement'
    TASKTITLE = 'Zanuda'
    TASKNAME = 'zanuda'
    TASKVERSION = 0.1
    DESCRIPTION = 'Space group validation'
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['aimless_pipe', 'coot_rebuild']

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData')

        self.createLine(['subtitle', 'Reflection data'])
        self.openSubFrame(frame=True)
        self.createLine(['widget', 'F_SIGF'])
        self.createLine(['widget', 'FREERFLAG'])
        self.closeSubFrame()

        self.createLine(['subtitle', 'Starting model'])
        self.openSubFrame(frame=True)
        self.createLine(['widget', 'XYZIN'])
        self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()

        self.createLine(['subtitle', 'Options'])
        self.openSubFrame(frame=True)
        self.createLine(
            [
                'widget',
                'AVERAGE',
                'label',
                'Symmetryse input model before further transformations'
            ]
        )
        self.closeSubFrame()
