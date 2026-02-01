import os
import pathlib

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class coordinate_selector(CPluginScript):
    TASKMODULE = 'model_data_utility'
    TASKTITLE = 'Select subset from coordinates'
    TASKNAME = 'coordinate_selector'
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'
    ERROR_CODES = { 202 : {'description' : 'Failed applying selection to PDB file' }
                    }
    PERFORMANCECLASS ='CAtomCountPerformance'

    def startProcess(self):
        try:
            self.container.inputData.XYZIN.loadFile()
            oldFullPath = pathlib.Path(self.container.outputData.XYZOUT.fullPath.__str__())
            if self.container.inputData.XYZIN.isMMCIF():
                self.container.outputData.XYZOUT.setFullPath(str(oldFullPath.with_suffix('.cif')))
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(str(self.container.outputData.XYZOUT.fullPath))
        except:
            raise
            self.appendErrorReport(202)
            return(CPluginScript.FAILED)
            
    def postProcessCheck(self, processId=None):
        if not os.path.isfile(str(self.container.outputData.XYZOUT.fullPath)):
            status = CPluginScript.FAILED
        else:
            status = CPluginScript.SUCCEEDED

        # Return tuple (status, exit_status, exit_code) as expected by base class
        exit_status = 0 if status == CPluginScript.SUCCEEDED else 1
        exit_code = 0  # No external process for this plugin
        return status, exit_status, exit_code
        
    def processOutputFiles(self):
        import gemmi

        if self.container.controlParameters.OVERRIDE_SUBTYPE.isSet():
            self.container.outputData.XYZOUT.subType = int(self.container.controlParameters.OVERRIDE_SUBTYPE)
        # Get the selection text, defaulting to "All" if not set or empty
        selection = self.container.inputData.XYZIN.selection
        if selection and selection.text and selection.text.isSet():
            selection_text = str(selection.text).strip()
            if not selection_text:
                selection_text = "All"
        else:
            selection_text = "All"
        self.container.outputData.XYZOUT.annotation.set(selection_text + ' of ' + self.container.inputData.XYZIN.annotation.__str__())

        from ccp4i2.core.CCP4ModelData import CPdbData
        aCPdbData = CPdbData()
        aCPdbData.loadFile(self.container.outputData.XYZOUT.fullPath)
        from lxml import etree
        rxml = etree.Element('CoordinateSelector')
        modelCompositionNode = etree.SubElement(rxml,'ModelComposition')

        try:
            st = gemmi.read_structure(self.container.outputData.XYZOUT.fullPath.__str__())
        except ValueError:
            try:
                st = gemmi.read_pdb(self.container.outputData.XYZOUT.fullPath.__str__())
            except:
                raise
        st.setup_entities()

        for mod in st:
            modEle = etree.SubElement(modelCompositionNode,"Model",id=str(mod.num))

        for model in st:
            for chain in model:
                for residue in chain:
                    if residue.entity_type == gemmi.EntityType.NonPolymer:
                        if len(residue) == 1 and residue[0].element.is_metal:
                            lig = etree.SubElement(modelCompositionNode,"Metal",id="/"+str(model.num)+"/"+chain.name+"/"+str(residue.seqid)+"("+residue.name+")")
                        else:
                            lig = etree.SubElement(modelCompositionNode,"Ligand",id="/"+str(model.num)+"/"+chain.name+"/"+str(residue.seqid)+"("+residue.name+")")

        model = st[0]

        nAtoms = 0
        nResidues = 0

        for chain in model:
            ch = etree.SubElement(modelCompositionNode,"Chain",id=chain.name)
            chain_name = etree.SubElement(ch,"Name")
            polylen = etree.SubElement(ch,"PolymerLength")
            waterlen = etree.SubElement(ch,"WaterLength")
            nonpolylen = etree.SubElement(ch,"NonPolymerLength")
            metallen = etree.SubElement(ch,"MetalLength")
            ligandlen = etree.SubElement(ch,"LigandLength")
            unknownlen = etree.SubElement(ch,"UnknownLength")
            polytype = etree.SubElement(ch,"PolymerTypes")
            chain_name.text = chain.name
            polylen.text = str(chain.get_polymer().length())
            waterlen.text = str(sum(res.is_water() for res in chain))
            nonpolylen.text = str(sum(res.entity_type == gemmi.EntityType.NonPolymer for res in chain))
            metallen.text = str(sum((res.entity_type == gemmi.EntityType.NonPolymer and len(res) == 1 and res[0].element.is_metal) for res in chain))
            ligandlen.text = str(sum((res.entity_type == gemmi.EntityType.NonPolymer and (len(res) != 1 or not res[0].element.is_metal)) for res in chain))
            unknownlen.text = str(sum(res.entity_type == gemmi.EntityType.Unknown for res in chain))
            polytypes = []
            for residue in chain:
                if residue.entity_type == gemmi.EntityType.Polymer:
                    nResidues += 1
                    nAtoms += len(residue)
                subchain = residue.subchain
                for entity in st.entities:
                    if subchain in entity.subchains:
                        if entity.polymer_type == gemmi.PolymerType.PeptideL:
                            polytypes.append("L-Peptide")
                        elif entity.polymer_type == gemmi.PolymerType.PeptideD:
                            polytypes.append("D-Peptide")
                        elif entity.polymer_type == gemmi.PolymerType.Pna:
                            polytypes.append("PNA")
                        elif entity.polymer_type == gemmi.PolymerType.Dna:
                            polytypes.append("DNA")
                        elif entity.polymer_type == gemmi.PolymerType.Rna:
                            polytypes.append("RNA")
                        elif entity.polymer_type == gemmi.PolymerType.DnaRnaHybrid:
                            polytypes.append("DNA/RNA")
                        elif entity.polymer_type == gemmi.PolymerType.SaccharideL:
                            polytypes.append("L-saccharide")
                        elif entity.polymer_type == gemmi.PolymerType.SaccharideD:
                            polytypes.append("D-saccharide")
                        elif entity.polymer_type == gemmi.PolymerType.CyclicPseudoPeptide:
                            polytypes.append("Cyclic pseudo-peptide")
            polytype.text = ",".join(list(set(polytypes)))

        with open(self.makeFileName('PROGRAMXML'),'w') as outputFile:
            CCP4Utils.writeXML(outputFile,etree.tostring(rxml, pretty_print=True))

        print("########################################")
        print("Set NATOMS")
        print("########################################")
        self.container.outputData.NATOMS.nAtoms.set(nAtoms)
        self.container.outputData.NATOMS.nResidues.set(nResidues)

        return CPluginScript.SUCCEEDED
