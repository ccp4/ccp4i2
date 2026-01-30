import os
import platform
import re
import sys

from lxml import etree

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript

from . import atomMatching, cifToMolBlock


class acedrgNew(CPluginScript):
    TASKMODULE = 'wrappers'                               # Where this plugin will appear on the gui
    TASKTITLE = 'acedrgNew'     # A short title for gui menu
    DESCRIPTION = 'Create a ligand dictionary with Acedrg'
    TASKNAME = 'acedrgNew'                                  # Task name - should be same as class name
    TASKCOMMAND = 'acedrg'                                     # The command to run the executable
    if platform.system() == 'Windows': TASKCOMMAND = 'acedrg.bat'
    TASKVERSION= 0.0                                     # Version of this plugin
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    ERROR_CODES = {  200 : { 'description' : 'Failed to add item to mol list' },201 : { 'description' : 'Failed to setFullPath' }, 202 : { 'description' : 'Failed to dump XMML' }, 203 : { 'description' : 'Failed to make RDKit Mol from DICT' },}
    
    def __init__(self,*args,**kws):
        CPluginScript.__init__(self, *args,**kws)
        self.xmlroot = etree.Element('Acedrg')
        self.smileStrCode = None

    def validity(self):
        """Filter CSMILESString validation errors when mode is not SMILES."""
        from ccp4i2.core import CCP4ErrorHandling
        error = super(acedrgNew, self).validity()
        mode = str(self.container.inputData.MOLORSMILES) if self.container.inputData.MOLORSMILES.isSet() else ""
        if mode != 'SMILES':
            # Filter out SMILES validation errors when not using SMILES input
            filtered = CCP4ErrorHandling.CErrorReport()
            for err in error.getErrors():
                err_class = err.get('class', '')
                err_name = err.get('name', '')
                if err_class == 'CSMILESString' and 'SMILESIN' in err_name:
                    continue
                filtered.append(
                    err.get('class', ''),
                    err.get('code', 0),
                    err.get('details', ''),
                    err.get('name', ''),
                    err.get('severity', 0)
                )
            return filtered
        return error

    def processInputFiles(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        try_mmCIF = False
        self.originalMolFilePath = None
        if self.container.inputData.MOLIN.isSet():
            self.originalMolFilePath = os.path.normpath(self.container.inputData.MOLIN.__str__())
        elif self.container.inputData.MOL2IN.isSet():
            self.originalMolFilePath = os.path.normpath(self.container.inputData.MOL2IN.__str__())

        elif self.container.inputData.MOLORSMILES.__str__() == 'SMILESFILE':
            self.originalMolFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'MOLIN.mol'))
            try:
                with open(self.container.inputData.SMILESFILEIN.__str__(),'r') as molinFile:
                    mol = molinFile.read().strip().split()[0]
                mol = Chem.MolFromSmiles(mol)
                AllChem.Compute2DCoords(mol)
                mol.SetProp("_MolFileChiralFlag","1")
                molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
                with open(self.originalMolFilePath,'w') as molinFile:
                    molinFile.write(molBlock)
            except:
                self.appendErrorReport(200,exc_info=sys.exc_info())
                return CPluginScript.FAILED

        elif self.container.inputData.MOLORSMILES.__str__() == 'SMILES':
            #Use rdkit to make mol from smiles...this will mean that we are always starting from a mol, and can
            #(hopefully) therefore assume that the atom order in our mol is the same as that in our dict
            self.originalMolFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'MOLIN.mol'))
            try:
                rdkinput = self.container.inputData.SMILESIN.__str__().partition(' ')[0]
                ssc_tmp = self.container.inputData.SMILESIN.__str__().partition(' ')[2].strip()
                if len(ssc_tmp) == 3:
                    self.smileStrCode = ssc_tmp
                mol = Chem.MolFromSmiles(rdkinput)
                AllChem.Compute2DCoords(mol)
                mol.SetProp("_MolFileChiralFlag","1")
                molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
                with open(self.originalMolFilePath,'w') as molinFile:
                    molinFile.write(molBlock)
            except:
                return CPluginScript.FAILED

        elif self.container.inputData.MOLORSMILES.__str__() == 'PDBMMCIF':
            import gemmi
            try:
                if gemmi.read_structure(self.container.inputData.PDBMMCIFIN.__str__()).input_format == gemmi.CoorFormat.Pdb:
                    self.originalMolFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'MOLIN.mol'))
                    mol = Chem.rdmolfiles.MolFromPDBFile(self.container.inputData.PDBMMCIFIN.__str__())
                    molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
                    with open(self.originalMolFilePath,'w') as molinFile:
                        molinFile.write(molBlock)
                else:
                    try_mmCIF = True
            except:
                self.appendErrorReport(200, exc_info=sys.exc_info())
                return CPluginScript.FAILED

        if self.container.inputData.MOLORSMILES.__str__() == 'DICT' or try_mmCIF == True:
            self.originalMolFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'MOLIN.mol'))
            print(self.originalMolFilePath)
            try:
                if self.container.inputData.DICTIN2.isSet():
                    molBlock = cifToMolBlock.cifFileToMolBlock(self.container.inputData.DICTIN2.__str__())
                elif self.container.inputData.PDBMMCIFIN.isSet() and try_mmCIF:
                    molBlock = cifToMolBlock.cifFileToMolBlock(self.container.inputData.PDBMMCIFIN.__str__())
                else:
                    pass #  should not happen
                print("molBlock:")
                print(molBlock)
                with open(self.originalMolFilePath,'w') as molinFile:
                    molinFile.write(molBlock)
            except:
                self.appendErrorReport(200, exc_info=sys.exc_info())
                return CPluginScript.FAILED
        #print 'Original mol file path', self.originalMolFilePath
    
        #Get the SMILES of the input MOL and put into report
        smilesNode = etree.SubElement(self.xmlroot,'SMILES')
        try:
            if self.container.inputData.MOL2IN.isSet():
                mol = Chem.MolFromMol2File(self.originalMolFilePath)
            else:
                mol = Chem.MolFromMolFile(self.originalMolFilePath)
            for iAtom in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(iAtom)
                atom.ClearProp('molAtomMapNumber')
            from rdkit.Chem import AllChem
            AllChem.Compute2DCoords(mol)
            smilesNode.text = Chem.MolToSmiles(mol, isomericSmiles=True)
        except:
            sys.stderr.write("ERROR: Getting SMILES from the input file failed.")
            smilesNode.text = "N/A"
        self.smilesString = smilesNode.text
            
        #I infer that ACEDRG reads in the bond and atom counts and the bonded atom indices as free input, whereas in practise
        #they are fixed format (i3,i3)...where nBonds > 99, this makes problems.  Insert a space to deal with this
        if os.path.isfile(self.originalMolFilePath):
            self.tmpMolFileName = os.path.normpath(os.path.join(self.getWorkDirectory(),'Kludged.mol'))
            with open(self.originalMolFilePath,'r') as molinFile:
                molBlock = molinFile.read()
                headerRead = False
                iAtom= 0
                iBond= 0
                molLines = molBlock.split('\n')
                for i in range(len(molLines)):
                    molLine = molLines[i]
                    if not headerRead:
                        if '999' in molLine.upper():
                            nBonds = int(molLine[3:6])
                            nAtoms = int(molLine[0:3])
                            if molLine[3:4] == ' ': newMolLine = molLine
                            else: newMolLine = molLine[0:3] + ' ' + molLine[3:]
                            molLines[i] = newMolLine
                            headerRead = True
                    elif iAtom < nAtoms:
                        iAtom += 1
                    elif iBond < nBonds:
                        if molLine[3:4] == ' ': newMolLine = molLine
                        else: newMolLine = molLine[0:3] + ' ' + molLine[3:]
                        molLines[i] = newMolLine
                        iBond += 1
                molBlock = '\n'.join(molLines)
                with open(self.tmpMolFileName,'w') as tmpMolFile:
                    tmpMolFile.write(molBlock)
    
    def makeCommandAndScript(self):
        if self.container.inputData.DICTIN2.isSet():
            self.appendCommandLine('-c')
            self.appendCommandLine(str(self.container.inputData.DICTIN2))
        elif self.container.inputData.MOLIN.isSet():
            self.appendCommandLine('-m')
            self.appendCommandLine(str(self.container.inputData.MOLIN))
        elif self.container.inputData.MOL2IN.isSet():
            self.appendCommandLine('-g')
            self.appendCommandLine(str(self.container.inputData.MOL2IN))
        elif self.container.inputData.PDBMMCIFIN.isSet():
            self.appendCommandLine('-x')
            self.appendCommandLine(str(self.container.inputData.PDBMMCIFIN))
        if self.container.inputData.METAL_STRUCTURE.isSet():
            self.appendCommandLine('--metalPDB=' + str(self.container.inputData.METAL_STRUCTURE))
        if self.container.controlParameters.NOPROT:
            self.appendCommandLine('--noProt')
        if self.container.inputData.DICTIN2.isSet() or self.container.inputData.MOLIN.isSet() or self.container.inputData.MOL2IN.isSet() or self.container.inputData.PDBMMCIFIN.isSet():
            if self.container.controlParameters.USE_COORD:
                self.appendCommandLine('-p')
        else:
            tmpSmileFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'tmp.smi'))
            with open(tmpSmileFilePath,'w') as tmpSmileFile:
                tmpSmileFile.write(self.smilesString)
            self.appendCommandLine('-i')
            self.appendCommandLine(tmpSmileFilePath)
        if self.container.inputData.NRANDOM and not self.container.controlParameters.USE_COORD:
            self.appendCommandLine('-j')
            self.appendCommandLine(self.container.inputData.NRANDOM)
        if self.container.inputData.TLC.__str__() or self.smileStrCode:
            self.appendCommandLine('-r')
            if self.smileStrCode != None and re.match("^[a-zA-Z0-9]*$", self.smileStrCode):
                self.appendCommandLine(self.smileStrCode)
            else:
                self.appendCommandLine(self.container.inputData.TLC.__str__())
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        print("acedrg is finished!!!")
        #Make database entry for Acedrg-generated PDB file
        pdbFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut.pdb'))
        self.container.outputData.XYZOUT.setFullPath(pdbFilePath)
        if self.smileStrCode != None and re.match("^[a-zA-Z0-9]*$", self.smileStrCode):
            annoCde = self.smileStrCode
        else:
            annoCde = self.container.inputData.TLC.__str__()
        retMatches = {}
        if os.path.isfile(pdbFilePath):
            try:
                import ccp4srs
                dummy = ccp4srs.Graph()
                print("Atom name matching available ...")
                if self.container.inputData.ATOMMATCHOPTION.__str__() == 'MONLIBCODE':
                    matchPdbFilePath = pdbFilePath[:-4]+"-matched.pdb"
                    retMatches = atomMatching.matchAtoms(pdbFilePath,matchPdbFilePath,self.container.inputData.MATCHTLC.__str__())
                    if len(retMatches) > 0:
                        self.container.outputData.XYZOUT.setFullPath(matchPdbFilePath)
                    print("Written",self.container.outputData.XYZOUT.fullPath.__str__())
                elif self.container.inputData.ATOMMATCHOPTION.__str__() == 'ALLMONLIB':
                    matchPdbFilePath = pdbFilePath[:-4]+"-matched.pdb"
                    retMatches = atomMatching.matchAtoms(pdbFilePath,matchPdbFilePath)
                    if len(retMatches) > 0:
                        self.container.outputData.XYZOUT.setFullPath(matchPdbFilePath)
                    print("Written",self.container.outputData.XYZOUT.fullPath.__str__())
                elif self.container.inputData.ATOMMATCHOPTION.__str__() == 'LOCALDICT':
                    matchPdbFilePath = pdbFilePath[:-4]+"-matched.pdb"
                    print("Supposed to be using local dictionary!!!!! Cannot yet!!!!!!",self.container.inputData.DICTIN.__str__())
                    retMatches = atomMatching.matchAtoms(pdbFilePath,matchPdbFilePath,dictFileName=self.container.inputData.DICTIN.__str__())
                    if len(retMatches) > 0:
                        self.container.outputData.XYZOUT.setFullPath(matchPdbFilePath)
                    print("Written",self.container.outputData.XYZOUT.fullPath.__str__())
                else:
                    self.container.outputData.XYZOUT.annotation = 'Coordinates for ligand '+annoCde
                    print("Written",self.container.outputData.XYZOUT.fullPath.__str__())
            except:
                #Atom naming requires more recent version of ccp4srs-python
                self.container.outputData.XYZOUT.annotation = 'Coordinates for ligand '+annoCde
                print("Written",self.container.outputData.XYZOUT.fullPath.__str__())

        #Make database entry for Acedrg-generated CIF file
        cifFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut.cif'))
        cifFilePath_alternative = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut_final.cif'))
        if not os.path.isfile(cifFilePath) and os.path.isfile(cifFilePath_alternative):
            cifFilePath = cifFilePath_alternative
        
        #Annoying kludge....for a brief window in time, Acedrg is not honouring the specified TLC
        if True:
            cifFilePathTmp = cifFilePath
            cifFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut_patched.cif'))
            with open(cifFilePathTmp,'r') as unpatchedFile:
                with open(cifFilePath,'w') as patchedFile:
                    lines = unpatchedFile.readlines()
                    for line in lines:
                        patchedFile.write(line.replace('UNL',annoCde))

        if os.path.isfile(cifFilePath):
            #Rename atoms in dictionary if we renamed them in PDB.
            if self.container.inputData.ATOMMATCHOPTION.__str__() != 'NOMATCH' and len(retMatches)>0:
                print("Matches !!!!!!!!!!",retMatches)
                print("Dict to match",cifFilePath)
                cifFilePathMatched = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut_matched_atoms.cif'))
                ls = atomMatching.replaceMatchesInDict(retMatches,cifFilePath,cifFilePathMatched)
                cifFilePath = cifFilePathMatched

            self.container.outputData.DICTOUT.setFullPath(cifFilePath)
            self.container.outputData.DICTOUT.annotation = 'Dictionary for ligand '+annoCde
            print("Written",self.container.outputData.DICTOUT.fullPath.__str__())
            
        #Here code to extract cif data into XML for representation as table
        if os.path.isfile(cifFilePath):
            try:
                from .MyCIFDigestor import MyCIFFile
                myCIFFile = MyCIFFile(filePath=cifFilePath)
                geometryNode = etree.SubElement(self.xmlroot,'Geometry')
                propertiesOfCategories = {'_chem_comp_bond':['atom_id_1','atom_id_2','value_dist','value_dist_esd','type'],
                '_chem_comp_angle':['atom_id_1','atom_id_2','atom_id_3','value_angle','value_angle_esd'],
                '_chem_comp_tor':['atom_id_1','atom_id_2','atom_id_3','atom_id_4','value_angle','period','value_angle_esd'],
                '_chem_comp_plane_atom':['atom_id','plane_id'],
                '_chem_comp_chir':['atom_id_centre','atom_id_1', 'atom_id_2', 'atom_id_3', 'volume_sign'  ],}
                blocks = [block for block in myCIFFile.blocks if block.category == 'data_comp_'+annoCde]
                for block in blocks:
                    for category in propertiesOfCategories:
                        properties = propertiesOfCategories[category]
                        examples = []
                        loops = [aLoop for aLoop in block.loops() if aLoop.category == category]
                        for loop in loops:
                            for iExample, loopline in enumerate(loop):
                                examples.append({'Number':str(iExample)})
                                for property in properties:
                                    if property in loopline:
                                        examples[-1][property] = loopline[property]
                        for example in examples:
                            exampleNode = etree.SubElement(geometryNode, category)
                            for property in properties:
                                if property in example:
                                    propertyNode = etree.SubElement(exampleNode, property)
                                    propertyNode.text = example[property]
            except:
                self.appendErrorReport(202)
                return CPluginScript.FAILED
            
        # AND THIS IS WHERE I START TRASHING STUFF ......

        from rdkit import Chem

        # Generate another RDKIT mol directly from the mol or smiles: this one *will* hopefully have proper
        # chirality information
        referenceMol = None
        referenceMolToDraw = None
        if self.container.inputData.MOL2IN.isSet():
            referenceMol = Chem.MolFromMol2File(self.originalMolFilePath)
            referenceMolToDraw = Chem.MolFromMol2File(self.originalMolFilePath)
        else:
            referenceMol = Chem.MolFromMolFile(self.originalMolFilePath)
            referenceMolToDraw = Chem.MolFromMolFile(self.originalMolFilePath)

        try:
            Chem.SanitizeMol(referenceMol)
            Chem.Kekulize(referenceMol)
            molToWrite = referenceMol
            # Output a MOL file
            molBlock = Chem.MolToMolBlock(molToWrite)
            with open(self.container.outputData.MOLOUT.fullPath.__str__(),'w') as outputMOL:
                outputMOL.write(molBlock)
                self.container.outputData.MOLOUT.annotation = 'MOL file from RDKIT'
            print("-----------------------")
            print("Written",self.container.outputData.MOLOUT.fullPath.__str__())
        except:
            import shutil
            molToWrite = referenceMol
            exc_type, exc_value = sys.exc_info()[:2]
            print("Failed writing MOL file")
            print(exc_type, exc_value)
            print(molToWrite)
            try:
                shutil.copy2(
                    os.path.normpath(os.path.join(self.getWorkDirectory(),'MOLIN.mol')),
                    self.container.outputData.MOLOUT.fullPath.__str__())
                print("Copied", os.path.normpath(os.path.join(self.getWorkDirectory(),'MOLIN.mol')),
                      "to", self.container.outputData.MOLOUT.fullPath.__str__())
            except:
                print("Failed copying MOL file")
            
        
        # Get 2D picture of structure from the RDKit mol and place in report
        svgNode = etree.SubElement(self.xmlroot,'SVGNode')
        try:
            from . import mol2svg
            svgText = bytes(mol2svg.svgFromMol(referenceMolToDraw),"utf-8")
            svgMolNode = etree.fromstring(svgText)
        except Exception as e:
            try:
                from ccp4i2.wrappers.Lidia.script import MOLSVG
                mdlMolecule = MOLSVG.MDLMolecule(self.container.outputData.MOLOUT.fullPath.__str__())
                svgMolNode = mdlMolecule.svgXML(size=(300,300))
            except Exception as e:
                print("ERROR: Drawing SVG picture of molecule was not successful.")
                print(e)
                svgMolNode = etree.fromstring("<svg></svg>")
        svgNode.append(svgMolNode)

        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot,pretty_print=True))

        return CPluginScript.SUCCEEDED
