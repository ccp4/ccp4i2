from ....qtgui import CCP4TaskWidget


class CTaskbuster(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'buster'
    TASKVERSION = 1.0
    TASKMODULE =['refinement']
    TASKTITLE = 'Refinement - BUSTER'
    SHORTTASKTITLE = "BUSTER"
    DESCRIPTION = 'Refine (BUSTER, Global Phasing Limited) with optional solvent/water update'
    PROGRAMHELP = ['buster']
    WHATNEXT = []

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData', title='Input Data')
        self.openSubFrame(frame=[True])
        self.createLine(['tip', 'Input reflections', 'widget', 'F_SIGF'])
        self.createLine(['tip', 'FreeR Flags', 'widget', 'FREERFLAG'])
        self.createLine(['tip', 'Model File', 'widget', 'XYZIN'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine(['tip', 'Ligand dictionary', 'widget', 'DICT'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine(['label', 'Perform', 'widget', 'NBCYCLES', 'label', 'refinements ("big cycles") with a maximum of', 'widget', 'NSCYCLES','label','iterations ("small cycles")'])
        self.closeSubFrame()
        self.setMenuText('WAT', {
                                 'ON' : 'Switch on water updates for every big cycle',
                                 'MAN': 'Switch on water solvent treatment for big cycle n onwards',
                                 'OFF': 'No automatic water update' })
        self.openSubFrame(frame=[True])
        self.createLine(['advice', 'Select Water Treatment' ])
        #self.createLine(['label', '   ', 'widget', '-guiMode', 'multiLineRadio', 'WAT'])
        self.createLine(['tip','BUSTER will use this option to determine how solvent is treated during refinement',
                         'label', 'Water solvent treatment:', 'widget', 'WAT'])  # 'stretch','widget'
        self.createLine(['label', 'From big cycle  ','widget','WATCYC','label', '  onwards'], toggle=['WAT', 'open', ['MAN']])
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine(['advice', 'Often used BUSTER options' ])
        self.createLine(['widget', 'AUTO_NCS', 'label', 'Automatically detect NCS (non-crystallographic symmetry) and impose LSSR (local structure similarity restraints)'])
        self.createLine(['widget', 'RBR', 'label', 'Turn on rigid body refinement for the first big cycle'])
        self.createLine(['widget', 'TLS', 'label', 'Carry out TLS (Translation/Libration/Screw) refinement using a simple model (one domain per chain)'])
        self.closeSubFrame()
