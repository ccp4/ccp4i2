"""
Martin Noble
"""

import os
import re

from PySide2 import QtCore
import gemmi

from ....core.CCP4Modules import PROJECTSMANAGER
from ....qtgui.CCP4TaskWidget import CTaskWidget


class CTaskProvideSequence(CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'ProvideSequence'
    TASKVERSION = 0.0
    TASKMODULE='data_entry'
    TASKTITLE="Import sequence(s)"
    WHATNEXT = []
    DESCRIPTION = '''Enter one or more sequences from a sequence file, from a PDB, or by cut and paste'''
    
    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)
        self.becauseSeqSet = False
        self.becauseAlignSet = False
#FIXME - SJM 2/6/2021 - I do not understand what this is for and it makes cooord/seqres choice tricky. I'm disabling it until I inderstand.
        self.becausePDBSet = False
 
    def drawContents(self):
        
        self.setProgramHelpFile('ProvideSequence')
        
        folder = self.openFolder(folderFunction='inputData',title='Optional objects from which to start definition')
        
        self.createLine( [ 'subtitle','Starting from file/object...'] )

        self.openSubFrame(frame=[True])
        self.createLine( [ 'widget', '-browseDb', True, '-enableEdit', False, 'SEQIN', 'tip', 'Sequence object or file' ] )
        self.container.inputData.SEQIN.dataChanged.connect( self.handleSelectSeqin)
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine( [ 'widget','-guiMode', 'radio', '-browseDb', True, 'XYZMODE', 'tip', 'Whether to determine sequence from coordinates or SEQRES' ])
        self.container.inputData.XYZMODE.dataChanged.connect(self.handleSelectXyzin)
        self.createLine( [ 'widget', '-browseDb', True, 'XYZIN', 'tip', 'Coordinate object or file' ] )
        self.container.inputData.XYZIN.dataChanged.connect( self.handleSelectXyzin)
        self.closeSubFrame()
        
        self.createLine( [ 'subtitle','Or enter the text of the sequence in a valid format...'] )
        self.createLine( [ 'widget', '-guiMode','multiLine','SEQUENCETEXT' ] )
        placeHolderText = """>IDENTIFIER A short description
PASTERYOURSEQUENCEINHERE"""
        self.getWidget('SEQUENCETEXT').widget.setPlaceholderText(placeHolderText)
        self.getWidget('SEQUENCETEXT').widget.textChanged.connect(self.validate)
#I am changing this 'unSet' to occur only if job is 'Pending'. This might still be unnecessary but better, I think.
        try:
            jobId = self.jobId()
            status = PROJECTSMANAGER().db().getJobInfo(jobId,'status')
            if status == "Pending":
                self.container.inputData.SEQIN.unSet()
        except:
            pass

    @QtCore.Slot()
    def handleSelectSeqin(self):
        #Here check to see whether this is a GUI-driven (rather than a programmatic) setting
        #the "self.becauseSeqSet" flag is set only when this is programmatic
        if self.becauseSeqSet: self.becauseSeqSet = False
        elif self.container.inputData.SEQIN.isSet():
            self.becauseSeqSet = True
            if os.path.isfile(self.container.inputData.SEQIN.fullPath.__str__()):
                with open(self.container.inputData.SEQIN.fullPath.__str__(),'r') as myFile:
                    content = myFile.read()
                    self.container.controlParameters.SEQUENCETEXT = content
            self.container.inputData.XYZIN.unSet()
            self.getWidget('SEQUENCETEXT').updateViewFromModel()
            self.getWidget('XYZIN').updateViewFromModel()
        self.validate()

    def sequencesFromSEQRES(self,fn):
        seqs = {}
        st = gemmi.read_structure(fn)
        st.setup_entities()
        for entity in st.entities:
            if entity.polymer_type == gemmi.PolymerType.PeptideL or entity.polymer_type == gemmi.PolymerType.PeptideD or entity.polymer_type == gemmi.PolymerType.Pna or entity.polymer_type == gemmi.PolymerType.Dna or entity.polymer_type == gemmi.PolymerType.Rna or entity.polymer_type == gemmi.PolymerType.DnaRnaHybrid or entity.polymer_type == gemmi.PolymerType.CyclicPseudoPeptide:
                seqs[str(entity.name)] = str(gemmi.one_letter_code(entity.full_sequence)+"\n")
        return seqs

    @QtCore.Slot()
    def handleSelectXyzin(self):
        if self.container.inputData.XYZIN.isSet():
            self.becausePDBSet = True
            if str(self.container.inputData.XYZMODE) == "seqres":
                sequences  = self.sequencesFromSEQRES( self.container.inputData.XYZIN.fullPath.__str__() )
            else:
                sequences  = self.sequencesFromPDB( self.container.inputData.XYZIN.fullPath.__str__() )
            self.container.controlParameters.SEQUENCETEXT = ''
            for sequenceId in sequences:
                if self.container.inputData.XYZIN.annotation.isSet():
                  formattedHeader = '>'+ re.sub(' ','_',str(self.container.inputData.XYZIN.annotation))
                else:
                  formattedHeader = '>'+os.path.split(self.container.inputData.XYZIN.fullPath.__str__())[1]
                formattedHeader+='_'+sequenceId+' Chain'+sequenceId + '\n'
                sequence = sequences[sequenceId]
                self.container.controlParameters.SEQUENCETEXT += formattedHeader
                self.container.controlParameters.SEQUENCETEXT += sequence
            self.container.inputData.SEQIN.unSet()
            self.getWidget('SEQUENCETEXT').updateViewFromModel()
            self.getWidget('SEQIN').updateViewFromModel()
        self.validate()

    def sequencesFromPDB(self, filePath):
        tlcOlcMap = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','TPO':'T','TYP':'Y','YPO':'Y'}
        waterNames = ['HOH','SOL','WAT','H2O']

        #FIXME - Should 5MV, OMG, etc be N rather than root NA?

        frequentSolutes = ['EDO','GOL','SUL','SO4','PO4','CL','ATP','PHO']
        tlcOlcMapNuc = { 'A':'A', 'ADE':'A', 'DA':'A', 'Ad':'A',
        'C':'C', 'CYT':'C', 'DC':'C', 'Cd':'C', 'DOC':'C', '5MC':'C', '5CM':'C', 'OMC':'C', '4OC':'C',
        'G':'G', 'GUA':'G', 'DG':'G', 'Gd':'G', '7MG':'G', 'OMG':'G', '2MG':'G', '8OG':'G', 
        'T':'T', 'THY':'T', 'DT':'T', 'Td':'T',
        'U':'U', 'URA':'U', '5MU': 'U', 'UR3':'U', 'H2U':'U', '4SU':'U', 'OMU':'U',
        'I':'I', 'INO':'I',
        'PSU':'Q',
        '3DR':'N', 'BRU':'N', 'MA6':'N',
        'DHU':'D',
        }
        sequences = {}
        if os.path.isfile(filePath):
            aCPdbData = CPdbData()
            aCPdbData.loadFile(filePath)
            mmdbManager = aCPdbData.mmdbManager
            for peptideChainId in aCPdbData.composition.peptides:
                sequences[peptideChainId] = ''
                peptideChain = mmdbManager.GetChain(1,peptideChainId)
                iOut = 0
                for i in range(peptideChain.GetNumberOfResidues()):
                    residueName = peptideChain.GetResidue(i).GetResName()
                    if residueName in tlcOlcMap:
                        sequences[peptideChainId]+= tlcOlcMap[residueName]
                        iOut += 1
                    elif (residueName in waterNames) or (residueName in frequentSolutes):
                        pass
                    else:
                        pass#print(residueName)
                        iOut += 1
                    if (iOut+1)%60 == 0 and sequences[peptideChainId][-1] != '\n':
                        sequences[peptideChainId]+= '\n'
                if sequences[peptideChainId][-1] != '\n': sequences[peptideChainId] += '\n'
            for nucleicChainId in aCPdbData.composition.nucleics:
                sequences[nucleicChainId] = ''
                nucleicChain = mmdbManager.GetChain(1,nucleicChainId)
                iOut = 0
                for i in range(nucleicChain.GetNumberOfResidues()):
                    residueName = nucleicChain.GetResidue(i).GetResName()
                    if residueName in tlcOlcMapNuc:
                        sequences[nucleicChainId]+= tlcOlcMapNuc[residueName]
                        iOut += 1
                    elif (residueName in waterNames) or (residueName in frequentSolutes):
                        pass
                    else:
                        print(residueName)
                        iOut += 1
                    if (iOut+1)%60 == 0 and sequences[nucleicChainId][-1] != '\n':
                        sequences[nucleicChainId]+= '\n'
                if sequences[nucleicChainId][-1] != '\n': sequences[nucleicChainId] += '\n'
                
        return sequences
    
    def isValid(self):
        #Here override logic of whether this is a valid task
        invalidElements = super(CTaskProvideSequence,self).isValid()
        if self.container.controlParameters.SEQUENCETEXT.isSet():
            if self.container.inputData.SEQIN in invalidElements:
                invalidElements.remove(self.container.inputData.SEQIN)
            if self.container.inputData.XYZIN in invalidElements:
                invalidElements.remove(self.container.inputData.XYZIN)
        return invalidElements
