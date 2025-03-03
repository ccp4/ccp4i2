from __future__ import print_function

import os
import sys
import networkx
from lxml import etree

import gemmi
import iotbx.phil
from iotbx.data_manager import DataManager
from mmtbx import process_predicted_model
from mmtbx.domains_from_pae import parse_pae_file

from core.CCP4PluginScript import CPluginScript

class editbfac(CPluginScript):

    TASKTITLE='Process Predicted Models'
    TASKNAME = 'editbfac'
    TASKMODULE=['alpha_fold', 'model_data_utility' ]
    TASKCOMMAND = 'None'
    TASKVERSION = 0.0
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS=False
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    def makeCommandAndScript(self):
        return

    def startProcess(self, comList, **kw):
        # Run cctbx conversion in startProcess. Setup cctbx dm & redirect std for this ftn.
        self.dm = DataManager()
        self.dm.set_overwrite(True)
        self.stdoutOrig = sys.stdout
        logfile = self.makeFileName('LOG') # from CCP4PluginScript
        sys.stdout = open(logfile, "w")
        # Setup parameters for cctbx (& PAE matrix from json file, if there is one)
        self.setupParams()
        self.setupPAE()
        # Prepare input & Load dist model (if there is one)
        inDistMod = self.container.inputData.XYZDISTMOD.fullPath.__str__()
        if os.path.isfile(inDistMod):
            self.distmod = self.dm.get_model(inDistMod)
        else:
            self.distmod = None
        print("========= AlphaFold/RosettaFold : i2 Process is converting pdb")
        runFile = self.setupInputModels()
        self.filelist = None
        self.convertFile(runFile)
        print("========= AF-RF Conversion done")
        sys.stdout.close()
        sys.stdout = self.stdoutOrig
    
    def setupInputModels(self):
        # Check & fix input files. Assumptions - a cif file will be an Alphafold 2 file (true when written). 
        # Robetta files can contain multiple models & will break the cctbx code (in 3/22).
        inFile = self.container.inputData.XYZIN.fullPath.__str__()
        cffin = gemmi.read_structure(inFile)
        # I assume the best choice is highest likelihood (standard convention). Remove rest.
        num_models = len(cffin)
        fparts = os.path.splitext(inFile)
        # Trouble. Robetta files are irregular (cctbx will reject them) & contain multiple models.
        # Currently can't access AUTHOR with Gemmi, & the REMARK's are stripped off ... which is unfortunate.
        isCifFile = os.path.splitext(os.path.split(inFile)[1])[1] == ".cif"
        isTRobFile = self.container.controlParameters.BTREATMENT.__str__() == "rmsd"
        # Fix the mess (keep the sep. in case I loop over the robetta models in the future).
        if num_models > 1:
            del cffin[1:num_models]
        if not (isTRobFile or isCifFile):
            if num_models > 1:
                nname = fparts[0] + "_alpconv.pdb"
                print("---> ALPHA PDB", nname)
                cffin.write_pdb(nname)
                return nname
            else:
                print("---> ALPHA PDB", inFile)
                return inFile
        if isCifFile:
            nname = fparts[0] + "_cifconv.pdb"
            print("---> ALPHA CIF", nname)
            cffin.write_pdb(nname)
            return nname
        if isTRobFile:
            nname = fparts[0] + "_robconv.pdb"
            print("---> ROSETTA", nname)
            cffin.write_pdb(nname)
            return nname

    def setupPAE(self):
        pae_file = self.container.inputData.PAEIN.fullPath.__str__()
        self.pae_matrix = None
        if os.path.isfile(pae_file):
            try:
                gopae = True
            except:
                gopae = False
                print("WARNING : networkx is not available in ccp4-python. Unable to process PAE Matrix in cctbx")
                print("You can install locally on Linux (Ubuntu) with ccp4-python -m pip install networkx")
            if gopae:
                try:
                    self.pae_matrix = parse_pae_file(pae_file)
                except:
                    print("WARNING : CCTBX failed to interpret the PAE file provided.")
                    print("Will proceed without PAE file.")

    def setupParams(self):
        master_phil = iotbx.phil.parse(process_predicted_model.master_phil_str)
        self.params = master_phil.extract()
        p = self.params.process_predicted_model
        # standard options
        p.b_value_field_is = self.container.controlParameters.BTREATMENT.__str__()  # 'plddt'
        p.remove_low_confidence_residues = self.container.controlParameters.CONFCUT      # True
        p.split_model_by_compact_regions = self.container.controlParameters.COMPACTREG   # True
        p.maximum_domains = self.container.controlParameters.MAXDOM  # 3
        p.domain_size = self.container.controlParameters.DOMAINSIZE.__float__()      # 15 Angstroms is this a float or what ? Not clear
        p.minimum_domain_length = self.container.controlParameters.MINDOML  # 10 nb. this is a float in the phil
        p.maximum_fraction_close = self.container.controlParameters.MAXFRACCL
        p.minimum_sequential_residues = self.container.controlParameters.MINSEQRESI
        p.minimum_remainder_sequence_length = self.container.controlParameters.MINREMSEQL
        p.minimum_plddt = self.container.controlParameters.MINLDDT.__float__() # 0.7
        p.maximum_rmsd = self.container.controlParameters.MAXRMSD.__float__() # 1.5
        # pae options
        p.pae_power = self.container.controlParameters.PAEPOWER
        p.pae_cutoff = self.container.controlParameters.PAECUTOFF
        p.pae_graph_resolution = self.container.controlParameters.PAEGRAPHRES
        # distance model options
        p.weight_by_ca_ca_distance = self.container.controlParameters.WEIGHTCA  # False
        p.distance_power = self.container.controlParameters.DISTPOW  # 1.0

    def convertFile(self, inFile):
        self.filelist = []
        cmfile = self.dm.get_model(inFile)
        model_info = process_predicted_model.process_predicted_model(cmfile, self.params, pae_matrix=self.pae_matrix,
                                                                     distance_model=self.distmod, log=sys.stdout)
        mmm = model_info.model.as_map_model_manager()
        # Prepare output pdb files (post translation)
        output_file_name = "converted_model.pdb"
        # print("CURRENTLY in :", os.getcwd(), " with file ", output_file_name, "  i2 CWD is ", self.getWorkDirectory())
        fofn = os.path.join( self.getWorkDirectory(), output_file_name)
        mmm.write_model(fofn)
        self.filelist.append(fofn)
        print("Writing model to file:- ", os.path.split(fofn)[1] )
        chainid_list = model_info.chainid_list
        if len(chainid_list) > 0:
            print("Model Segments found: %s" %(" ".join(chainid_list)))
            for chainid in chainid_list:
                selection_string = "chain %s" %(chainid)
                ph = model_info.model.get_hierarchy()
                asc1 = ph.atom_selection_cache()
                sel = asc1.selection(selection_string)
                m1 = model_info.model.select(sel)
                outp = os.path.join( self.getWorkDirectory(), '%s_chain%s.pdb' %(output_file_name[:-4], chainid))
                print("Writing chain %s to file:- "%(chainid), os.path.split(outp)[1])
                self.filelist.append(outp)
                self.dm.write_model_file(m1, outp, chainid)

    def processOutputFiles(self):
        status = CPluginScript.FAILED
        outputXYZFILES = self.container.outputData.XYZFILES
        for afile in self.filelist:
            outputXYZFILES.append(outputXYZFILES.makeItem())
            outputXYZFILES[-1].setFullPath(afile)
            outputXYZFILES[-1].annotation = os.path.basename(afile)
        status = CPluginScript.SUCCEEDED
        # Create a trivial xml output file
        from core import CCP4File
        root = etree.Element('editbfac')
        self.container.outputData.XYZOUT.subType = 1
        f = CCP4File.CXmlDataFile(fullPath=self.makeFileName('PROGRAMXML'))
        f.saveFile(root)
        return status

# Stuart's previous stuff - changed.
#===============================================================================
#     def convert_plddt_to_bfactor(self, struct):
#         """
#         Adam's code from MrParse to convert pLDDT scores to bfactors based on the original cctbx function (process_predicted_model)
#         """
#         for chain in struct[0]:
#             for residue in chain:
#                 for atom in residue:
#                     plddt_value = atom.b_iso
#                     atom.b_iso = self._convert_to_bfactor(plddt=plddt_value)
#         return struct
#     
#     def convert_rmsd_to_bfactor(self, struct):
#         """
#         Convert rmsd scores to bfactors based on the original cctbx function (process_predicted_model)
#         """
#         for chain in struct[0]:
#             for residue in chain:
#                 for atom in residue:
#                     rmsd_value = atom.b_iso
#                     atom.b_iso = self._convert_to_bfactor(rmsd=rmsd_value)
#         return struct
#     
#     def _convert_to_bfactor(self, plddt=None, rmsd=None):
#         """
#         pLDDT/rmsd conversion algorithm
#         """
#         if plddt is not None:
#             lddt = plddt / 100
#             if lddt <= 0.5:
#                 return 657.97  # Same as the b-factor value with an rmsd estimate of 5.0
#             rmsd_est = (0.6 / (lddt ** 3))
#         elif plddt is None and rmsd is not None:
#             # Make sure error estimates are in reasonable range
#             if rmsd < 0:
#                 rmsd = 0
#             elif rmsd > 20:
#                 rmsd = 20
#             rmsd_est = rmsd
#         elif plddt is None and rmsd is None:
#             rmsd_est = 1.2
#         bfactor = ((8 * (np.pi ** 2)) / 3.0) * (rmsd_est ** 2)
#         if bfactor > 999.99:
#             return 999.99
#         return bfactor
# 
#     def alphaFoldAutoSet(self,st):
#         """
#         Call to Adam's convertor for pLDDT values
#         """
#         self.convert_plddt_to_bfactor(st)
# 
#     def rosettaFoldAutoSet(self,st):
#         """
#         Call to Adam's convertor for rmsd values
#         """
#         self.convert_rmsd_to_bfactor(st)
# 
#     def startProcess(self,comList,**kw):
#     
#         inp = self.container.inputData
#         par = self.container.controlParameters
#         out = self.container.outputData
# 
#         import gemmi
# 
#         inputFile = inp.XYZIN.fullPath.__str__()
#         outputFile = out.XYZOUT.fullPath.__str__()
# 
#         st = gemmi.read_structure(inputFile)
# 
#         if par.BTREATMENT.isSet() :
#             if str( par.BTREATMENT ) == 'a' :
#                 self.alphaFoldAutoSet(st)
#             elif str( par.BTREATMENT ) == 'r' :
#                self.rosettaFoldAutoSet(st)
#             elif str( par.BTREATMENT ) == 'f' :
#                 for model in st:
#                     for cra in model.all():
#                         cra.atom.b_iso = par.BFACTOR_STATIC
#             else :
#                raise Exception()
#         else:
#            raise Exception()
# 
# #        if par.AUTOSETALPHAFOLD:
# #            self.alphaFoldAutoSet(st)
# #        elif par.AUTOSETROSETTAFOLD:
# #            self.rosettaFoldAutoSet(st)
# #        else:
# #            for model in st:
# #                for cra in model.all():
# #                    cra.atom.b_iso = par.BFACTOR_STATIC
# 
#         st.write_pdb(outputFile)
# 
#     def processOutputFiles(self):
#         import os
#         status = CPluginScript.FAILED
#         if os.path.exists(self.container.outputData.XYZOUT.__str__()): status = CPluginScript.SUCCEEDED
#         print('editbfac.handleFinish',self.container.outputData.XYZOUT, os.path.exists(self.container.outputData.XYZOUT.__str__()))
#         # Create a trivial xml output file
#         from core import CCP4File
#         from lxml import etree
#         root = etree.Element('editbfac')
# 
#         self.container.outputData.XYZOUT.subType = 1
# 
#         f = CCP4File.CXmlDataFile(fullPath=self.makeFileName('PROGRAMXML'))
#         f.saveFile(root)
#         return status
#===============================================================================

