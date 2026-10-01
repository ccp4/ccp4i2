import os
import sys

import gemmi
from lxml import etree
try:  # iotbx/mmtbx (cctbx) present only in the execution (worker) env, not the slim API
    import iotbx.phil
    from iotbx.data_manager import DataManager
    from mmtbx import process_predicted_model
    from mmtbx.domains_from_pae import parse_pae_file
except ImportError:
    iotbx = DataManager = process_predicted_model = parse_pae_file = None

from ccp4i2.core.CCP4PluginScript import CPluginScript


def residue_ranges(numbers):
    """[26, 27, 28, 40] -> "26-28, 40"."""
    ranges, start, prev = [], None, None
    for n in sorted(numbers):
        if start is None:
            start = prev = n
        elif n == prev + 1:
            prev = n
        else:
            ranges.append((start, prev))
            start = prev = n
    if start is not None:
        ranges.append((start, prev))
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in ranges)


def describe_model(path):
    """(number of residues, their ranges) of a model file's first model."""
    structure = gemmi.read_structure(str(path))
    numbers = [r.seqid.num for c in structure[0] for r in c] if len(structure) else []
    return len(numbers), residue_ranges(numbers)


def summarise(input_path, model_path, domain_paths):
    """What was kept, as an XML element: the input's residue count, the
    processed model's, and each domain's chain, residues and ranges."""
    root = etree.Element("editbfac")
    n_in, _ = describe_model(input_path)
    etree.SubElement(root, "InputResidues").text = str(n_in)
    n_kept, kept = describe_model(model_path)
    model = etree.SubElement(root, "Model", residues=str(n_kept), ranges=kept)
    model.text = os.path.basename(model_path)
    for path in domain_paths:
        n, ranges = describe_model(path)
        chain = os.path.splitext(os.path.basename(path))[0].rsplit("_chain", 1)[-1]
        domain = etree.SubElement(root, "Domain", chain=chain, residues=str(n), ranges=ranges)
        domain.text = os.path.basename(path)
    return root


class editbfac(CPluginScript):
    TASKNAME = 'editbfac'

    def startProcess(self):
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
        # Plain Python values, not the parameters themselves. (They were passed
        # as CFloat/CInt/CBoolean objects: 1 / pae ** CFloat raised inside
        # cctbx's PAE clustering, which swallows the exception and returns
        # None, so every job given a PAE file failed with "'NoneType' object
        # is not iterable".)
        c = self.container.controlParameters
        val = lambda item, kind: kind(item.value)
        # standard options
        p.b_value_field_is = str(c.BTREATMENT)  # 'plddt'
        p.remove_low_confidence_residues = val(c.CONFCUT, bool)
        p.split_model_by_compact_regions = val(c.COMPACTREG, bool)
        p.maximum_domains = val(c.MAXDOM, int)
        p.domain_size = val(c.DOMAINSIZE, float)
        p.minimum_domain_length = val(c.MINDOML, float)
        p.maximum_fraction_close = val(c.MAXFRACCL, float)
        p.minimum_sequential_residues = val(c.MINSEQRESI, int)
        p.minimum_remainder_sequence_length = val(c.MINREMSEQL, int)
        p.minimum_plddt = val(c.MINLDDT, float)
        p.maximum_rmsd = val(c.MAXRMSD, float)
        # pae options
        p.pae_power = val(c.PAEPOWER, float)
        p.pae_cutoff = val(c.PAECUTOFF, float)
        p.pae_graph_resolution = val(c.PAEGRAPHRES, float)
        # distance model options
        p.weight_by_ca_ca_distance = val(c.WEIGHTCA, bool)
        p.distance_power = val(c.DISTPOW, float)

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
        # Say what each file holds (they were listed by file name alone, so a
        # user could not tell which domain was which without opening them),
        # and record it for the report.
        summary = summarise(self.container.inputData.XYZIN.fullPath.__str__(),
                            self.filelist[0], self.filelist[1:])
        n_in = summary.findtext("InputResidues")
        model = summary.find("Model")
        domains = {d.text: d for d in summary.findall("Domain")}
        outputXYZFILES = self.container.outputData.XYZFILES
        for afile in self.filelist:
            outputXYZFILES.append(outputXYZFILES.makeItem())
            outputXYZFILES[-1].setFullPath(afile)
            d = domains.get(os.path.basename(afile))
            if d is None:
                text = "Processed model: %s of %s residues kept" % (model.get("residues"), n_in)
            else:
                text = "Domain %s: residues %s (%s residues)" % (
                    d.get("chain"), d.get("ranges"), d.get("residues"))
            outputXYZFILES[-1].annotation = text
        from ccp4i2.core import CCP4File
        self.container.outputData.XYZOUT.subType = 1
        f = CCP4File.CXmlDataFile(fullPath=self.makeFileName('PROGRAMXML'))
        f.saveFile(summary)
        return CPluginScript.SUCCEEDED
