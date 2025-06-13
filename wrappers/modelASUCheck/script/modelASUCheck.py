"""
    modelASUCheck.py: CCP4 GUI Project

    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.

    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

import os
import traceback
from lxml import etree
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils

def addSequenceAlignmentToEtree(fn,asuin,root):
    import gemmi
    import numpy as np
    from scipy.optimize import linear_sum_assignment

    xml_seq_align = etree.SubElement(root,"SequenceAlignment")

    provide_seq_asu = []

    for idx in range(len(asuin.fileContent.seqList)):
        seq = asuin.fileContent.seqList[idx].sequence
        nCopies = int(asuin.fileContent.seqList[idx].nCopies)
        for ic in range(nCopies):
            if asuin.fileContent.seqList[idx].polymerType == "PROTEIN":
                seq_full = [gemmi.expand_one_letter(i, gemmi.ResidueKind.AA) for i in seq]
                provide_seq_asu.append((seq_full,gemmi.ResidueKind.AA))
            elif asuin.fileContent.seqList[idx].polymerType == "DNA":
                seq_full = [gemmi.expand_one_letter(i, gemmi.ResidueKind.DNA) for i in seq]
                provide_seq_asu.append((seq_full,gemmi.ResidueKind.DNA))
            elif asuin.fileContent.seqList[idx].polymerType == "RNA":
                seq_full = [gemmi.expand_one_letter(i, gemmi.ResidueKind.RNA) for i in seq]
                provide_seq_asu.append((seq_full,gemmi.ResidueKind.RNA))

    st = gemmi.read_structure(fn)
    st.setup_entities()

    align_scores = []
    align_results = []
    for c in st[0]:
        score_row = []
        results_row = []
        for asu_seq,seq_type in provide_seq_asu:
            if seq_type == gemmi.ResidueKind.AA:
                matchType = gemmi.PolymerType.PeptideL
            elif seq_type == gemmi.ResidueKind.DNA:
                matchType = gemmi.PolymerType.Dna
            elif seq_type == gemmi.ResidueKind.RNA:
                matchType = gemmi.PolymerType.Rna
            result = gemmi.align_sequence_to_polymer(asu_seq,
                             c.get_polymer(),
                             matchType,
                             gemmi.AlignmentScoring())

            score_row.append(result.score)
            results_row.append((result,asu_seq,c))

        align_scores.append(score_row)
        align_results.append(results_row)

    cost = np.array(align_scores)
    row_ind, col_ind = linear_sum_assignment(cost,True)

    ii = 0
    for icol in col_ind:
        irow = row_ind[ii]
        result = align_results[irow][icol][0]
        seq =  align_results[irow][icol][1]
        chain =  align_results[irow][icol][2]
        ii += 1

        xml_this_chain_align = etree.SubElement(xml_seq_align,"Alignment")
        xml_chain_id = etree.SubElement(xml_this_chain_align,"ChainID")
        xml_chain_id.text = chain.name
        xml_chain_match_count = etree.SubElement(xml_this_chain_align,"match_count")
        xml_chain_match_count.text = str(result.match_count)
        xml_chain_match_identity = etree.SubElement(xml_this_chain_align,"identity")
        xml_chain_match_identity.text = str(result.calculate_identity())
        xml_chain_match_identity_1 = etree.SubElement(xml_this_chain_align,"identity_1")
        xml_chain_match_identity_1.text = str(result.calculate_identity(1))
        xml_chain_match_identity_2 = etree.SubElement(xml_this_chain_align,"identity_2")
        xml_chain_match_identity_2.text = str(result.calculate_identity(2))
        xml_chain_match_cigar = etree.SubElement(xml_this_chain_align,"CIGAR")
        xml_chain_match_cigar.text = str(result.cigar_str())

        xml_chain_match_align_1 = etree.SubElement(xml_this_chain_align,"align_1")
        xml_chain_match_align_1.text = result.add_gaps(gemmi.one_letter_code(seq),1)
        xml_chain_match_align_1 = etree.SubElement(xml_this_chain_align,"align_match")
        xml_chain_match_align_1.text = result.match_string
        xml_chain_match_align_2 = etree.SubElement(xml_this_chain_align,"align_2")
        xml_chain_match_align_2.text = result.add_gaps(gemmi.one_letter_code(chain.get_polymer().extract_sequence()), 2)

class modelASUCheck(CPluginScript):
    TASKNAME = 'modelASUCheck'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'person@server.com'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="None"
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS = False

    def __init__(self, *args, **kws):
        super(modelASUCheck, self).__init__(*args, **kws)

    def makeCommandAndScript(self,**kw):
        return
    
    def startProcess(self, comList, **kw):
    
        xyzin = self.container.inputData.XYZIN
        asuin = self.container.inputData.ASUIN

        self.xmlroot = etree.Element('modelASUCheck')
        
        if asuin.isSet() and xyzin.isSet():
            try:
                addSequenceAlignmentToEtree(str(xyzin),asuin,self.xmlroot)
            except Exception as err:
                traceback.print_exc()
                print("...importing sequences for alignment test failed", err)
        CCP4Utils.saveEtreeToFile(self.xmlroot,self.makeFileName('PROGRAMXML'))

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        return CPluginScript.SUCCEEDED
