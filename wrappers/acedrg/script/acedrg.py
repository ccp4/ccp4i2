from __future__ import print_function


from core.CCP4PluginScript import CPluginScript
from PySide6 import QtCore
import os,glob,re,time,sys
from core import CCP4XtalData
from lxml import etree
import math
from core import CCP4Modules,CCP4Utils
import platform

class acedrg(CPluginScript):
    TASKMODULE = 'wrappers'                               # Where this plugin will appear on the gui
    TASKTITLE = 'acedrg'     # A short title for gui menu
    DESCRIPTION = 'Sketch a ligand'
    TASKNAME = 'acedrg'                                  # Task name - should be same as class name
    TASKCOMMAND = 'acedrg'                                     # The command to run the executable
    if platform.system() == 'Windows': TASKCOMMAND = 'acedrg.bat'
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    #WHATNEXT = ['coot_rebuild','parrot','buccaneer_build_refine_mr']
    #RUNEXTERNALPROCESS=False

    ERROR_CODES = {  200 : { 'description' : 'Failed to add item to mol list' },201 : { 'description' : 'Failed to setFullPath' }, 202 : { 'description' : 'Failed to dump XMML' }, 203 : { 'description' : 'Failed to make RDKit Mol from DICT' },}
    
    def __init__(self,*args,**kws):
        CPluginScript.__init__(self, *args,**kws)
        self.xmlroot = etree.Element('Acedrg')
    
    def processInputFiles(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        if self.container.inputData.MOLIN.isSet():
            self.originalMolFilePath = os.path.normpath(self.container.inputData.MOLIN.__str__())
        if self.container.inputData.MOLORSMILES.__str__() == 'SMILES':
            #Use rdkit to make mol from smiles...this will mean that we are always starting from a mol, and can
            #(hopefully) therefore assume that the atom order in our mol is the same as that in our dict
            self.originalMolFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'MOLIN.mol'))
            try:
                mol = Chem.MolFromSmiles(self.container.inputData.SMILESIN.__str__())
                AllChem.Compute2DCoords(mol)
                mol.SetProp("_MolFileChiralFlag","1")
                molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
                with open(self.originalMolFilePath,'w') as molinFile:
                    molinFile.write(molBlock)
            except:
                return CPluginScript.FAILED
        #print 'Original mol file path', self.originalMolFilePath
    
        #Get the SMILES of the input MOL and put into report
        mol = Chem.MolFromMolFile(self.originalMolFilePath)
        for iAtom in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(iAtom)
            atom.ClearProp('molAtomMapNumber')
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
        smilesNode = etree.SubElement(self.xmlroot,'SMILES')
        smilesNode.text = Chem.MolToSmiles(mol, isomericSmiles=True)
        self.smilesString = smilesNode.text
        
        #I infer that ACEDRG reads in the bond and atom counts and the bonded atom indices as free input, whereas in practise
        #they are fixed format (i3,i3)...where nBonds > 99, this makes problems.  Insert a space to deal with this
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
        tmpSmileFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'tmp.smi'))
        with open(tmpSmileFilePath,'w') as tmpSmileFile:
            tmpSmileFile.write(self.smilesString)
        #self.appendCommandLine('-m')
        #self.appendCommandLine(self.tmpMolFileName)
        self.appendCommandLine('-z')
        self.appendCommandLine('-i')
        self.appendCommandLine(tmpSmileFilePath)
        self.appendCommandLine('-r')
        self.appendCommandLine(self.container.inputData.TLC.__str__())
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        #Make database entry for Acedrg-generated PDB file
        pdbFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut.pdb'))
        if os.path.isfile(pdbFilePath):
            self.container.outputData.XYZOUT.setFullPath(pdbFilePath)
            self.container.outputData.XYZOUT.annotation = 'Coordinates for ligand '+self.container.inputData.TLC.__str__()

        #Make database entry for Acedrg-generated CIF file
        cifFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut.cif'))
        
        #Annuying kludge....for a brief window in time, Acedrg is not honouring the specified TLC
        if True:
            cifFilePathTmp = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut.cif'))
            cifFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'AcedrgOut_patched.cif'))
            with open(cifFilePathTmp,'r') as unpatchedFile:
                with open(cifFilePath,'w') as patchedFile:
                    lines = unpatchedFile.readlines()
                    for line in lines:
                        patchedFile.write(line.replace('UNL',self.container.inputData.TLC.__str__()))
            
        if os.path.isfile(cifFilePath):
            self.container.outputData.DICTOUT.setFullPath(cifFilePath)
            self.container.outputData.DICTOUT.annotation = 'Dictionary for ligand '+self.container.inputData.TLC.__str__()
            
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
                blocks = [block for block in myCIFFile.blocks if block.category == 'data_comp_'+self.container.inputData.TLC.__str__()]
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
                self.apppendErrorReport(202)
                return CPluginScript.FAILED
            
        from rdkit import Chem
        # Generate RDKit mol from the Dict: this one won't necessarily (or even probably) have proper carry over
        # of chirality information
        mol = molFromDict(cifFilePath)
        molToWrite = mol

        # Generate another RDKIT mol directly from the mol or smiles: this one *will* hopefully have proper
        # chirality information
        referenceMol = None
        referenceMol = Chem.MolFromMolFile(self.originalMolFilePath)
        Chem.SanitizeMol(referenceMol)
        Chem.Kekulize(referenceMol)

        #Find largest common substructure, and corresponding atom equivalences
        from rdkit.Chem import MCS
        mols = [mol, referenceMol]
        atomsOfMol = [i for i in range(mol.GetNumAtoms())]
        atomsOfTemplate = [i for i in range(mol.GetNumAtoms())]
        try:
            MCSResult = MCS.FindMCS(mols,bondCompare='any')
            pattern = Chem.MolFromSmarts(MCSResult.smarts)
            atomsOfMol = mol.GetSubstructMatch(pattern)
            atomsOfTemplate = referenceMol.GetSubstructMatch(pattern)
        except:
            warningNode = etree.SubElement(self.xmlroot,'Warning')
            warningNode.text = 'RDKIT failed to use MaximumCommonStructure matching - chirality less likely to be right'
        for iEquiv in range(len(atomsOfMol)):
            atomOfMol = mol.GetAtomWithIdx(atomsOfMol[iEquiv])
            atomOfTemplate = referenceMol.GetAtomWithIdx(atomsOfTemplate[iEquiv])
            atomOfTemplate.SetMonomerInfo(atomOfMol.GetMonomerInfo())
        molToWrite = referenceMol


        # Generate a low energy PDB file from the RDKit mol
        pdbBlock, sortedEnergyArray = lowEnergyPDBForMol(molToWrite,int(self.container.inputData.NRANDOM))
        with open(self.container.outputData.XYZOUT_RDKIT.fullPath.__str__(),'w') as outputRDKitPDB:
            outputRDKitPDB.write(pdbBlock)
            self.container.outputData.XYZOUT_RDKIT.annotation = 'Coordinates from RDKIT'
        iConf = 0
        
        #generate a mol from the pdb block
        from rdkit.Chem.rdmolfiles import MolFromPDBBlock
        pdbMol = MolFromPDBBlock(pdbBlock)
    
        #Provide data to allow graphing of conformer energies
        for cid, energy in sortedEnergyArray:
            confNode = etree.SubElement(self.xmlroot,'Conformer')
            rankNode = etree.SubElement(confNode,'Rank')
            rankNode.text = str(iConf)
            energyNode = etree.SubElement(confNode,'Energy')
            energyNode.text = str(energy)
            iConf += 1
        
        # Output a MOL file
        molBlock = Chem.MolToMolBlock(molToWrite)
        with open(self.container.outputData.MOLOUT.fullPath.__str__(),'w') as outputMOL:
            outputMOL.write(molBlock)
            self.container.outputData.MOLOUT.annotation = 'MOL file from RDKIT'

        #Here code to make a "patched" CIF file with atom coordinates from the molToWrite
        if os.path.isfile(cifFilePath):
            if True:
                from .MyCIFDigestor import MyCIFFile
                myCIFFile = MyCIFFile(filePath=cifFilePath)
                blocks = [block for block in myCIFFile.blocks if block.category == 'data_comp_'+self.container.inputData.TLC.__str__()]
                bestConformer = pdbMol.GetConformer(id=0)
                for block in blocks:
                    loops = [aLoop for aLoop in block.loops() if aLoop.category == '_chem_comp_atom']
                    for loop in loops:
                        for iAtom in range(pdbMol.GetNumAtoms()):
                            atom = pdbMol.GetAtomWithIdx(iAtom)
                            atomName = atom.GetMonomerInfo().GetName()
                            atomPosition = bestConformer.GetAtomPosition(iAtom)
                            #Here identify loop lines that correspond to atom with this name
                            looplines = [aLoopline for aLoopline in loop if aLoopline['atom_id'] == atomName.strip()]
                            for loopline in looplines:
                                print('A match ',atomName.strip(), ' ',loopline['atom_id'])
                                loop.setLinePropertyValue(loopline,'x',str(atomPosition.x))
                                #Here I update the value in loopline, since this is used as a key in setLinePropertyValue, so must reflect contents thereof
                                loopline['x'] = str(atomPosition.x)
                                loop.setLinePropertyValue(loopline,'y',str(atomPosition.y))
                                loopline['y'] = str(atomPosition.y)
                                loop.setLinePropertyValue(loopline,'z',str(atomPosition.z))
                                loopline['z'] = str(atomPosition.z)
                        #Now set all hydrogen atom positions to be unknown
                        looplines = [aLoopline for aLoopline in loop if aLoopline['type_symbol'] == 'H']
                        for loopline in looplines:
                            loop.setLinePropertyValue(loopline,'x','.')
                            #Here I update the value in loopline, since this is used as a key in setLinePropertyValue, so must reflect contents thereof
                            loopline['x'] = '.'
                            loop.setLinePropertyValue(loopline,'y','.')
                            loopline['y'] = '.'
                            loop.setLinePropertyValue(loopline,'z','.')
                            loopline['z'] = '.'
                with open (self.container.outputData.DICTOUT_RDKIT.fullPath.__str__(),'w') as dictoutRdkit:
                    dictoutRdkit.write(str(myCIFFile))
                self.container.outputData.DICTOUT_RDKIT.annotation = 'Acedrg restraints, RDKit conformer for '+self.container.inputData.TLC.__str__()
            else:
                pass

        # Get 2D picture of structure from the RDKit mol and place in report
        svgNode = etree.SubElement(self.xmlroot,'SVGNode')
        svgNode.append(svgFromMol(molToWrite))

        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot, pretty_print=True))

        return CPluginScript.SUCCEEDED

def svgFromMol(mol):
    try:
        from rdkit.Chem.Draw import spingCanvas
        import rdkit.Chem.Draw.MolDrawing
        from rdkit.Chem.Draw.MolDrawing import DrawingOptions
        
        myCanvas = spingCanvas.Canvas(size=(350,350),name='MyCanvas',imageType='svg')
        myDrawing = rdkit.Chem.Draw.MolDrawing(canvas=myCanvas)
        for iAtom in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(iAtom)
            atom.ClearProp('molAtomMapNumber')
            '''
            print 'ANr',atom.GetAtomicNum()
            print 'FC',atom.GetFormalCharge()
            print 'NRadEl',atom.GetNumRadicalElectrons()
            print 'Isot',atom.GetIsotope()
            print 'HasProp',atom.HasProp('molAtomMapNumber')
            print 'Degree',atom.GetDegree()
        print 'noCarbSym',myDrawing.drawingOptions.noCarbonSymbols
        print 'includeAtom',myDrawing.drawingOptions.includeAtomNumbers'''
        myDrawing.AddMol(mol)
        svg = myCanvas.canvas.text()#.replace('svg:','')
        return etree.fromstring(svg)
    except:
        from rdkit import Chem
        molBlock = Chem.MolToMolBlock(mol)
        import MOLSVG
        mdlMolecule = MOLSVG.MDLMolecule(molBlock=molBlock)
        return mdlMolecule.svgXML(size=(350,350))

def svgFromPDB(pdbFilePath):
    from rdkit import Chem
    mol = Chem.rdmolfiles.MolFromPDBFile(pdbFilePath)
    Chem.Kekulize(mol)
    return svgFromRdkitMol(mol)
    
def svgFromRdkitMol(mol):
    
    from rdkit.Chem.Draw import spingCanvas
    import rdkit.Chem.Draw.MolDrawing
    from rdkit.Chem.Draw.MolDrawing import DrawingOptions
    
    myCanvas = spingCanvas.Canvas(size=(300,300),name='MyCanvas',imageType='svg')
    myDrawing = rdkit.Chem.Draw.MolDrawing(canvas=myCanvas)
    
    
    from rdkit.Chem.Draw.MolDrawing import DrawingOptions
    drawingOptions = DrawingOptions()
    drawingOptions.noCarbonSymbols=True
    drawingOptions.includeAtomNumbers=False
    myDrawing.AddMol(mol, drawingOptions=drawingOptions)
    svg = myCanvas.canvas.text().replace('svg:','')
    return etree.fromstring(svg)


def molFromDict(cifFilePath):
    from rdkit import Chem
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    mmcif_dict = MMCIF2Dict ( cifFilePath)

    from core import CCP4ModelData
    elMap = {}
    for iElement in range (len(CCP4ModelData.CElement.QUALIFIERS['enumerators'])):
        elMap[CCP4ModelData.CElement.QUALIFIERS['enumerators'][iElement].upper()] = iElement+1
    bondTypeMap = {'aromatic':Chem.rdchem.BondType.AROMATIC,'single':Chem.rdchem.BondType.SINGLE,'double':Chem.rdchem.BondType.DOUBLE,'triple':Chem.rdchem.BondType.TRIPLE,'1.5':Chem.rdchem.BondType.ONEANDAHALF,'deloc':Chem.rdchem.BondType.ONEANDAHALF}
    atomIdMap = {}

    eMol = Chem.EditableMol(Chem.Mol())

    conformer = Chem.Conformer(len(mmcif_dict['_chem_comp_atom.type_symbol']))

    for iAtom, type_symbol in enumerate(mmcif_dict['_chem_comp_atom.type_symbol']):
        type_symbol = type_symbol.upper()
        raw_atom_id = mmcif_dict['_chem_comp_atom.atom_id'][iAtom].strip()
        atom_id = ''
        if len(type_symbol) == 1: atom_id += ' '
        atom_id += raw_atom_id
        while len(atom_id)<4: atom_id += ' '
        comp_id = mmcif_dict['_chem_comp_atom.comp_id'][iAtom]
        atomPDBResidueInfo = Chem.rdchem.AtomPDBResidueInfo()
        atomPDBResidueInfo.SetResidueName(comp_id)
        atomPDBResidueInfo.SetResidueNumber(1)
        atomPDBResidueInfo.SetName(atom_id)
        atomPDBResidueInfo.SetOccupancy(1.0)
        atomPDBResidueInfo.SetTempFactor(20.0)
        atomPDBResidueInfo.SetIsHeteroAtom(True)
        atom = Chem.Atom(elMap[type_symbol])
        atom.SetProp('molAtomMapNumber',atom_id)
        atom.SetMonomerInfo(atomPDBResidueInfo)
        atomIdMap[raw_atom_id] = iAtom
        from rdkit.Geometry.rdGeometry import Point3D
        try:
            xString = mmcif_dict['_chem_comp_atom.x'][iAtom]
            yString = mmcif_dict['_chem_comp_atom.y'][iAtom]
            zString = mmcif_dict['_chem_comp_atom.z'][iAtom]
        except:
            try:
                xString = mmcif_dict['_chem_comp_atom.pdbx_model_Cartn_x_ideal'][iAtom]
                yString = mmcif_dict['_chem_comp_atom.pdbx_model_Cartn_y_ideal'][iAtom]
                zString = mmcif_dict['_chem_comp_atom.pdbx_model_Cartn_z_ideal'][iAtom]
            except:
                xString = '.'
                yString = '.'
                zString = '.'
        if xString.strip() != '.' and yString.strip() != '.' and zString.strip() != '.':
            atomPosition = Point3D(float(xString),
                                   float(yString),
                                   float(zString),)
            conformer.SetAtomPosition(iAtom,atomPosition)
        eMol.AddAtom(atom)

    for iBond in range(len(mmcif_dict['_chem_comp_bond.atom_id_2'])):
        atom1Number = atomIdMap[mmcif_dict['_chem_comp_bond.atom_id_1'][iBond]]
        atom2Number = atomIdMap[mmcif_dict['_chem_comp_bond.atom_id_2'][iBond]]
        try:
            bondType = bondTypeMap[mmcif_dict['_chem_comp_bond.type'][iBond].lower()]
        except KeyError as err:
            print('Could not match bond type', mmcif_dict['_chem_comp_bond.type'][iBond].lower())
            bondType = Chem.rdchem.BondType.SINGLE
        except:
            bondType = Chem.rdchem.BondType.SINGLE
            if mmcif_dict['_chem_comp_bond.value_order'][iBond] == 'DOUB': bondType = Chem.rdchem.BondType.DOUBLE
        #print('Adding bond', atom1Number, atom2Number, bondType)
        eMol.AddBond(atom1Number,atom2Number,bondType)

    mol = eMol.GetMol()
    initialConformerId = mol.AddConformer(conformer,assignId=True)
    blk = Chem.MolToPDBBlock(mol)
    mol = Chem.MolFromPDBBlock(blk,sanitize=True)

    from rdkit.Chem import Draw
    try:
        Chem.SanitizeMol(mol)
    except:
        print('Sanitize mol did not work')

    # Use the coordinates in the input CIF file to assign stereochemistry
    try:
        #print(initialConformerId)
        from rdkit.Chem import rdmolops
        rdmolops.AssignAtomChiralTagsFromStructure(mol, confId=initialConformerId, replaceExistingTags=True)
    except:
        print('Unable to assign stereochemistry from coordinates in the dict file')

    try:
        Chem.Kekulize(mol)
    except:
        pass

    from rdkit.Chem import AllChem
    confId2D = AllChem.Compute2DCoords(mol)
    return mol

def lowEnergyPDBForMol(molInput, nRandom=500):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(molInput)
    cids = AllChem.EmbedMultipleConfs(mol, nRandom, clearConfs = False)
    
    #Here something for clustering would be desirable
    conformerEnergyArray = []
    for icid, cid in enumerate(cids):
        if icid != 0:
            AllChem.UFFOptimizeMolecule(mol, confId=cid)
            forceField = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            confEnergy = forceField.CalcEnergy()
            pair = (cid, confEnergy)
            conformerEnergyArray.append(pair)
    sortedArray = sorted(conformerEnergyArray, key=lambda conformer: conformer[1])
    for icid, cid in enumerate(cids):
        if mol.GetConformer(id=cid).Is3D() and cid != sortedArray[0][0]:
            mol.RemoveConformer(cid)

    molRemovedHs = Chem.RemoveHs(mol)
    pdbBlock = Chem.MolToPDBBlock(molRemovedHs, sortedArray[0][0])

    return pdbBlock, sortedArray
