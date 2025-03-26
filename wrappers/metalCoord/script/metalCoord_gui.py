from qtgui import CCP4TaskWidget

class CmetalCoord_gui(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'metalCoord'
    TASKVERSION = 0.2
    TASKMODULE =['refinement']
    TASKTITLE = 'MetalCoord'
    SHORTTASKTITLE = 'MetalCoord'
    DESCRIPTION = 'Generate restraints for metal-containing monomers'
    WHATNEXT = ['servalcat_pipe']

    def __init__(self, parent):
        CCP4TaskWidget.CTaskWidget.__init__(self, parent)

    def drawContents(self):
        tab = '&nbsp;&nbsp;&nbsp;&nbsp;'
        self.openFolder(folderFunction='inputData', title='Input Data')
        self.createLine(['tip', 'Atomic model', 'widget', 'XYZIN'])
        self.createLine(['label', 'Monomer code', 'widget', 'LIGAND_CODE'])
        self.openSubFrame(frame=[True], title='Advanced parameters')
        self.createLine( [ 'label', 'Distance threshold: (range 0-1)<br/><i>A threshold d to select atoms is (r<sub>1</sub> + r<sub>2</sub>)*(1 + d) where r<sub>1</sub> and r<sub>2</sub> are covalent radii.</i>', 'stretch', 'widget', 'DISTANCE_THRESHOLD'])

        self.createLine( [ 'label', 'Maximum coordination number:', 'stretch', 'widget', 'MAXIMUM_COORDINATION_NUMBER'])
        for i in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,24]:
            self.createLine( [ 'label', 'Coordination class:', 'stretch', 'widget', f'COORD{i:02d}'], toggle=['MAXIMUM_COORDINATION_NUMBER', 'open', [str(i)]])

        self.createLine( [ 'label', 'Procrustes distance threshold: (range 0-1)', 'stretch', 'widget', 'PROCRUSTES_DISTANCE_THRESHOLD'])
        self.createLine( [ 'label', 'Minimum sample size for statistics:', 'stretch', 'widget', 'MINIMUM_SAMPLE_SIZE'])
        self.createLine( [ 'widget', 'USE_PDB', 'label', 'Use COD structures based on the input PDB/mmCIF coordinates'])
        self.createLine( [ 'widget', 'IDEAL_ANGLES', 'label', 'Provide only ideal bond angles'])
        self.createLine( [ 'widget', 'SIMPLE', 'label', 'Simple distance based filtering'])
        self.createLine( [ 'widget', 'SAVE_PDBMMCIF', 'label', 'Update link records to metal sites in the atomic model'])
        self.createLine( [ 'widget', 'KEEP_LINKS', 'label', 'Keep existing link records to metal sites'], toggle=['SAVE_PDBMMCIF', 'open', [True]])
        self.closeSubFrame()
