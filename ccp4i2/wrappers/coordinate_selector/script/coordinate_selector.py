import os
import pathlib
import xml.etree.ElementTree as ET

import gemmi

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class coordinate_selector(CPluginScript):
    
    TASKMODULE = 'model_data_utility' # Where this plugin will appear on gui
    TASKTITLE = 'Select subset from coordinates' # A short title for gui menu
    TASKNAME = 'coordinate_selector'   # Task name - should be same as class name
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection to PDB file' }
                    }
    TASKCOMMAND='None'
    RUNEXTERNALPROCESS=False
    PERFORMANCECLASS ='CAtomCountPerformance'
    
    def makeCommandAndScript(self,**kw):
      return
    
    def startProcess(self,comList,**kw):
        try:
            self.container.inputData.XYZIN.loadFile()
            oldFullPath = pathlib.Path(self.container.outputData.XYZOUT.fullPath.__str__())
            if self.container.inputData.XYZIN.isMMCIF():
                self.container.outputData.XYZOUT.setFullPath(str(oldFullPath.with_suffix('.cif')))
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(str(self.container.outputData.XYZOUT.fullPath))
        except:
            raise

    def postProcessCheck(self, processId):
        if not os.path.isfile(str(self.container.outputData.XYZOUT.fullPath)): return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
        
    def processOutputFiles(self):
        if self.container.controlParameters.OVERRIDE_SUBTYPE.isSet():
            self.container.outputData.XYZOUT.subType = int(self.container.controlParameters.OVERRIDE_SUBTYPE)
        self.container.outputData.XYZOUT.annotation.set(self.container.inputData.XYZIN.selection.__str__()+' of '+self.container.inputData.XYZIN.annotation.__str__())

        from ....core.CCP4ModelData import CPdbData
        aCPdbData = CPdbData()
        aCPdbData.loadFile(self.container.outputData.XYZOUT.fullPath)
        rxml = ET.Element('CoordinateSelector')
        modelCompositionNode = ET.SubElement(rxml,'ModelComposition')

        try:
            st = gemmi.read_structure(self.container.outputData.XYZOUT.fullPath.__str__())
        except ValueError:
            try:
                st = gemmi.read_pdb(self.container.outputData.XYZOUT.fullPath.__str__())
            except:
                raise
        st.setup_entities()

        for mod in st:
            modEle = ET.SubElement(modelCompositionNode,"Model",id=str(mod.num))

        for model in st:
            for chain in model:
                for residue in chain:
                    if residue.entity_type == gemmi.EntityType.NonPolymer:
                        if len(residue) == 1 and residue[0].element.is_metal:
                            lig = ET.SubElement(modelCompositionNode,"Metal",id="/"+str(model.num)+"/"+chain.name+"/"+str(residue.seqid)+"("+residue.name+")")
                        else:
                            lig = ET.SubElement(modelCompositionNode,"Ligand",id="/"+str(model.num)+"/"+chain.name+"/"+str(residue.seqid)+"("+residue.name+")")

        model = st[0]

        nAtoms = 0
        nResidues = 0

        for chain in model:
            ch = ET.SubElement(modelCompositionNode,"Chain",id=chain.name)
            chain_name = ET.SubElement(ch,"Name")
            polylen = ET.SubElement(ch,"PolymerLength")
            waterlen = ET.SubElement(ch,"WaterLength")
            nonpolylen = ET.SubElement(ch,"NonPolymerLength")
            metallen = ET.SubElement(ch,"MetalLength")
            ligandlen = ET.SubElement(ch,"LigandLength")
            unknownlen = ET.SubElement(ch,"UnknownLength")
            polytype = ET.SubElement(ch,"PolymerTypes")
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

        CCP4Utils.writeXml(rxml, self.makeFileName('PROGRAMXML'))

        print("########################################")
        print("Set NATOMS")
        print("########################################")
        self.container.outputData.NATOMS.nAtoms.set(nAtoms)
        self.container.outputData.NATOMS.nResidues.set(nResidues)

        return CPluginScript.SUCCEEDED
