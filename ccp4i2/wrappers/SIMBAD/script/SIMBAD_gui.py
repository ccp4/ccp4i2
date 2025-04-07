from multiprocessing import cpu_count

from ....qtgui.CCP4TaskWidget import CTaskWidget


class SIMBAD_gui(CTaskWidget):
    TASKNAME = 'SIMBAD' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = 'molecular_replacement'
    SHORTTASKTITLE='SIMBAD Molecular Replacement Pipeline'
    TASKTITLE='Sequence Free Molecular Replacement - SIMBAD'
    DESCRIPTION = 'This task is for running Molecular Replacement without a sequence'
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['prosmart_refmac', 'coot_rebuild']

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    def drawContents(self):

        self.setProgramHelpFile('SIMBAD')
        self.openFolder(folderFunction='inputData', followFrom=False)

        self.openSubFrame(toggle=['PERFORM', 'close', ['den']])
        self.createLine(['advice', 'Input reflections'])
        self.createLine(['widget', 'F_SIGF'])
        self.closeSubFrame()

        self.createLine(['subtitle', 'Search level:', 'widget', 'SIMBAD_SEARCH_LEVEL'])
        self.createLine(['subtitle', 'Organism:', 'widget', 'SIMBAD_ORGANISM'], toggle=['SIMBAD_SEARCH_LEVEL', 'close',

        # Number of processors
        self.openSubFrame(frame=True)
        # There's no sensible way to set dynamic defaults so we use an unusual number or processors in the xml to indicate we haven't set a dynamic default yet
        if self.container.inputData.SIMBAD_NPROC == 9993:  self.container.inputData.SIMBAD_NPROC = cpu_count()
        x = self.container.inputData.SIMBAD_NPROC.qualifiers()['guiLabel']
        self.createLine(['subtitle', x, 'widget', 'SIMBAD_NPROC'])
        self.closeSubFrame()
        
        self.drawOptions()
    
    def drawOptions(self):
        self.openFolder(folderFunction='inputData',title='Advanced Options')
        label = self.container.inputData.SIMBAD_ROT_PROGRAM.qualifiers()['guiLabel']
        self.createLine(['subtitle', label, 'widget', 'SIMBAD_ROT_PROGRAM'], toggle=['SIMBAD_SEARCH_LEVEL', 'close',
                                                                                     ['Lattice']])
        label = self.container.inputData.SIMBAD_MR_PROGRAM.qualifiers()['guiLabel']
        self.createLine(['subtitle', label, 'widget', 'SIMBAD_MR_PROGRAM'])
        label = self.container.inputData.SIMBAD_NMOL.qualifiers()['guiLabel']
        self.createLine(['subtitle', label, 'widget', 'SIMBAD_NMOL'])
        label = self.container.inputData.SIMBAD_PROCESS_ALL.qualifiers()['guiLabel']
        self.createLine(['subtitle', label, 'widget', 'SIMBAD_PROCESS_ALL'])
        label = self.container.inputData.SIMBAD_SGALTERNATIVE.qualifiers()['guiLabel']
        self.createLine(['subtitle', label, 'widget', 'SIMBAD_SGALTERNATIVE'])
