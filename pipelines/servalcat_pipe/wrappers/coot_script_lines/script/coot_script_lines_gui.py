"""
    tasks/coot_script_lines
    """

from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class Ccoot_script_lines(CTaskWidget):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'coot_script_lines'
    TASKVERSION = 0.1
    TASKMODULE='model_building'
    TASKTITLE  = 'Scripted model building - COOT'
    SHORTTASKTITLE  = 'Scripted COOT'
    DESCRIPTION='''Use scripts to fit sidechains, perform stepped refinement, fill and fit... (non-interactive Coot)'''
    
    def drawContents(self):
        
        
        self.setProgramHelpFile('coot_script_lines')
        
        self.openFolder(folderFunction='inputData')
        
        self.createLine( [ 'widget', 'XYZIN' ] )
        self.createLine( [ 'subtitle', 'Electron density maps:' ] )
        self.createLine( [ 'widget', 'FPHIIN' ] )
        self.createLine( [ 'subtitle', 'Difference density maps:' ] )
        self.createLine( [ 'widget', 'DELFPHIIN' ] )
        self.createLine( [ 'widget', 'DICT' ] )
        
        #self.openFolder(folderFunction='controlParameters',title='Options')
        self.createLine(['advice'])
        self.createLine(['subtitle','Coot operation to perform'])
        self.createLine( [ 'widget', 'STARTPOINT'] )
        self.container.controlParameters.STARTPOINT.dataChanged.connect(self.handleStartpointChanged)
        self.createLine(['label','Show/edit corresponding coot script','widget','SHOWSCRIPT'])
        self.createLine( [ 'widget', '-guiMode','multiLine','SCRIPT' ],toggle=['SHOWSCRIPT','open',[True]] )

    @QtCore.Slot()
    def handleStartpointChanged(self):
        if str(self.container.controlParameters.STARTPOINT) == 'BLANK':
            self.container.controlParameters.SHOWSCRIPT=True
            self.container.controlParameters.SCRIPT='''# Input molecules are given internal identifiers "MolHandle_1" etc
# Input maps are given internal identifiers "MapHandle_1" etc
# Input difference maps are given internal identifiers "DifmapHandle_1" etc
#
# The appropriate place to put output pdb files to be picked up by the system
# is stored in the variable dropDir (see below for how to use this)
#
# Beware spaces...this is python after all
#
# Example commands you may wish to try can be explored by uncommenting the lines below
# fill_partial_residues(MolHandle_1)
# fit_protein(MolHandle_1)
# stepped_refine_protein_for_rama(MolHandle_1)
#
# Per chain operations might look something like:
# chainAMol = new_molecule_by_atom_selection(MolHandle_1,'//A')
# otherChains = ['B','C','D','E','F','G','H','I','J','K','L']
# for otherChain in otherChains:
#    molToMove = new_molecule_by_atom_selection(MolHandle_1,'//A')
#    newMolHandle = superpose_with_chain_selection(MolHandle_1, molToMove, otherChain, 'A', 1, 1, 1)
#    merge_molecules(newMolHandle, molToMove)
#
# The following line ensures the output gets written out:
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif str(self.container.controlParameters.STARTPOINT) == 'FILL_PARTIAL_RESIDUES':
            self.container.controlParameters.SCRIPT='''# Input molecules are given internal identifiers "MolHandle_1" etc"
# Input maps are given internal identifiers "MapHandle_1" etc"
# Input difference maps are given internal identifiers "DifmapHandle_1" etc"
#
# The appropriate place to put output pdb files to be picked up by the system
# is stored in the variable dropDir (see below for how to use this)
#
# Beware spaces...this is python after all
fill_partial_residues(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif str(self.container.controlParameters.STARTPOINT) == 'FIT_PROTEIN':
            self.container.controlParameters.SCRIPT='''# Input molecules are given internal identifiers "MolHandle_1" etc"
# Input maps are given internal identifiers "MapHandle_1" etc"
# Input difference maps are given internal identifiers "DifmapHandle_1" etc"
#
# The appropriate place to put output pdb files to be picked up by the system
# is stored in the variable dropDir (see below for how to use this)
#
# Beware spaces...this is python after all
fit_protein(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif str(self.container.controlParameters.STARTPOINT) == 'STEPPED_REFINE_PROTEIN_FOR_RAMA':
            self.container.controlParameters.SCRIPT='''# Input molecules are given internal identifiers "MolHandle_1" etc"
# Input maps are given internal identifiers "MapHandle_1" etc"
# Input difference maps are given internal identifiers "DifmapHandle_1" etc"
#
# The appropriate place to put output pdb files to be picked up by the system
# is stored in the variable dropDir (see below for how to use this)
#
# Beware spaces...this is python after all
stepped_refine_protein_for_rama(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif str(self.container.controlParameters.STARTPOINT) == 'STEPPED_REFINE_PROTEIN':
            self.container.controlParameters.SCRIPT='''# Input molecules are given internal identifiers "MolHandle_1" etc"
# Input maps are given internal identifiers "MapHandle_1" etc"
# Input difference maps are given internal identifiers "DifmapHandle_1" etc"
#
# The appropriate place to put output pdb files to be picked up by the system
# is stored in the variable dropDir (see below for how to use this)
#
# Beware spaces...this is python after all
stepped_refine_protein(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif str(self.container.controlParameters.STARTPOINT) == 'MORPH_FIT':
            self.container.controlParameters.SCRIPT='''# Input molecules are given internal identifiers "MolHandle_1" etc"
# Input maps are given internal identifiers "MapHandle_1" etc"
# Input difference maps are given internal identifiers "DifmapHandle_1" etc"
#
# The appropriate place to put output pdb files to be picked up by the system
# is stored in the variable dropDir (see below for how to use this)
#
# Beware spaces...this is python after all
def is_polymer_chain(imol, ch_id):
    return is_protein_chain_qm(imol, ch_id) or is_nucleotide_chain_qm(imol, ch_id)


def reso_morph(imol, imol_map, n_rounds):
    for round in range(n_rounds):
        f = float(round)/float(n_rounds)
        for ch_id in chain_ids(imol):
            if is_polymer_chain(imol, ch_id):

                # play with these numbers
                radius =   6 * (2 - f)
                sf     = 100 * (1 - f)
                sharpen(0, sf)
                morph_fit_chain(imol, ch_id, radius)

# run the script
reso_morph(MolHandle_1, MapHandle_1, 20)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
