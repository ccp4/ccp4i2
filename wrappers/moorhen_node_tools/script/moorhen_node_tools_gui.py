"""
tasks/moorhen_node_tools
"""

from PySide2 import QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class Cmoorhen_node_tools(CTaskWidget):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'moorhen_node_tools'
    TASKVERSION = 0.1
    TASKMODULE='model_building'
    TASKTITLE  = 'Scripted model building - Moorhen Headless'
    SHORTTASKTITLE  = 'Scripted Moorhen Headless'
    DESCRIPTION='''Use scripts to fit sidechains, perform stepped refinement, fill and fit... (non-interactive Moorhen)'''
    
    def drawContents(self):
        
        
        self.setProgramHelpFile('moorhen_node_tools')
        
        folder = self.openFolder(folderFunction='controlParameters',title='Basic data')
        self.openSubFrame( frame=[True] )
        self.createLine(['subtitle','Coot operation to perform'])
        self.createLine( [ 'widget', 'STARTPOINT'] )
        self.closeSubFrame()
        
        self.openSubFrame( frame=[True] )
        self.createLine( [ 'subtitle', 'Models:' ] )
        self.createLine( [ 'widget', 'XYZIN' ] )
        self.createLine( [ 'subtitle', 'Density map coefficients:' ] )
        self.createLine( [ 'widget', 'FPHIIN' ] )
        self.createLine( [ 'subtitle', 'Difference density map coefficients:' ] )
        self.createLine( [ 'widget', 'DELFPHIIN' ] )
        self.createLine( [ 'subtitle', 'Density maps:' ] )
        self.createLine( [ 'widget', 'MAPIN' ] )
        self.createLine( [ 'subtitle', 'Restraint dictionaries:' ] )
        self.createLine( [ 'widget', 'DICTIN' ] )
        self.closeSubFrame()
        
    @QtCore.Slot()
    def handleStartpointChanged(self):
        if str(self.container.controlParameters.STARTPOINT) == 'BLANK':
            self.container.controlParameters.SHOWSCRIPT=True
            self.container.controlParameters.SCRIPT='''
import os
import coot_headless_api as chapi

mc = chapi.molecules_container_t(True)

def int_list_to_mh_string(solutions):
    s = ""
    for ii,imol in enumerate(solutions):
        if ii == 0:
            s += str(imol)
        else:
            s += ":" + str(imol)
    return s

def fit_ligs(MolHandle_1, MapHandle_1, TLC):

    imol_lig = mc.get_monomer(TLC)
    if mc.is_valid_model_molecule(MolHandle_1):
        if mc.is_valid_map_molecule(MapHandle_1):
            solutions = mc.fit_ligand(MolHandle_1, MapHandle_1, imol_lig, 1.0, True, 30)
            print("Found", len(solutions), "solutions")
            if solutions:
                sl = int_list_to_mh_string(solutions)
                mc.merge_molecules(MolHandle_1, sl)

def exec(plugin):
    mc.read_pdb(plugin.container.inputData.XYZIN[0].__str__())
    mc.read_mtz(plugin.container.inputData.FPHIIN[0].__str__(), 'F', 'PHI', '', False, False)
    
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
