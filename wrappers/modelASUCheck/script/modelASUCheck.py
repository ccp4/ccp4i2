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

import traceback

from lxml import etree
from scipy.optimize import linear_sum_assignment
import gemmi
import numpy as np

from core import CCP4Utils
from core.CCP4PluginScript import CPluginScript


class modelASUCheck(CPluginScript):
    TASKNAME = "modelASUCheck"
    TASKVERSION = 0.1
    MAINTAINER = "person@server.com"
    ERROR_CODES = {
        201: {"description": "Failed to analyse output files"},
        202: {"description": "Failed applying selection ot PDB file"},
    }
    RUNEXTERNALPROCESS = False

    def startProcess(self, command=None, handler=None, **kw):
        xyzin = self.container.inputData.XYZIN
        asuin = self.container.inputData.ASUIN

        xml = etree.Element("modelASUCheck")
        try:
            xml.append(sequenceAlignment(str(xyzin), asuin))
        except Exception as err:
            traceback.print_exc()
            print("...importing sequences for alignment test failed", err)
            return CPluginScript.FAILED
        CCP4Utils.saveEtreeToFile(xml, self.makeFileName("PROGRAMXML"))
        return CPluginScript.SUCCEEDED


def sequenceAlignment(xyzinPath, asuin):
    structure = gemmi.read_structure(xyzinPath)
    structure.setup_entities()

    align_scores = []
    align_results = []
    for chain in structure[0]:
        score_row = []
        results_row = []
        for fullSeq, polymerType in sequences(asuin):
            scoring = gemmi.AlignmentScoring()
            result = gemmi.align_sequence_to_polymer(
                fullSeq, chain.get_polymer(), polymerType, scoring
            )
            score_row.append(result.score)
            results_row.append((result, fullSeq, chain))
        align_scores.append(score_row)
        align_results.append(results_row)

    cost = np.array(align_scores)
    row_ind, col_ind = linear_sum_assignment(cost, True)

    matched_model_chains = []
    xml = etree.Element("SequenceAlignment")
    for irow, icol in zip(row_ind, col_ind):
        result, seq, chain = align_results[irow][icol]
        matched_model_chains.append(chain.name)
        alignment = etree.SubElement(xml, "Alignment")
        etree.SubElement(alignment, "ChainID").text = chain.name
        etree.SubElement(alignment, "match_count").text = str(result.match_count)
        etree.SubElement(alignment, "identity").text = str(result.calculate_identity())
        etree.SubElement(alignment, "identity_1").text = str(result.calculate_identity(1))
        etree.SubElement(alignment, "identity_2").text = str(result.calculate_identity(2))
        etree.SubElement(alignment, "CIGAR").text = str(result.cigar_str())
        etree.SubElement(alignment, "align_1").text = result.add_gaps(gemmi.one_letter_code(seq), 1)
        etree.SubElement(alignment, "align_match").text = result.match_string
        etree.SubElement(alignment, "align_2").text = result.add_gaps(
            gemmi.one_letter_code(chain.get_polymer().extract_sequence()), 2
        )

    polymer_types = [gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD, gemmi.PolymerType.Dna, gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid]

    nonAligned = etree.SubElement(xml, "NonAlignedModelChains")

    list_seq = []
    for seq in asuin.fileContent.seqList:
        residueKind, polymerType = {
            "PROTEIN": (gemmi.ResidueKind.AA, gemmi.PolymerType.PeptideL),
            "DNA": (gemmi.ResidueKind.DNA, gemmi.PolymerType.Dna),
            "RNA": (gemmi.ResidueKind.RNA, gemmi.PolymerType.Rna),
        }[seq.polymerType]
        list_seq.append((gemmi.expand_one_letter_sequence(str(seq.sequence), residueKind),polymerType))

    for chain in structure[0]:
        if not chain.name in matched_model_chains and chain.get_polymer().check_polymer_type() in polymer_types:

            iseq = 0
            for fullSeq, polymerType in list_seq:
                alignment = etree.SubElement(nonAligned, "Alignment")
                scoring = gemmi.AlignmentScoring()
                result = gemmi.align_sequence_to_polymer(fullSeq, chain.get_polymer(), polymerType, scoring)
                etree.SubElement(alignment, "SeqAUNo").text = str(iseq)
                etree.SubElement(alignment, "ChainID").text = chain.name
                etree.SubElement(alignment, "match_count").text = str(result.match_count)
                etree.SubElement(alignment, "identity").text = str(result.calculate_identity())
                etree.SubElement(alignment, "identity_1").text = str(result.calculate_identity(1))
                etree.SubElement(alignment, "identity_2").text = str(result.calculate_identity(2))
                etree.SubElement(alignment, "CIGAR").text = str(result.cigar_str())
                etree.SubElement(alignment, "align_1").text = result.add_gaps(gemmi.one_letter_code(fullSeq), 1)
                etree.SubElement(alignment, "align_match").text = result.match_string
                etree.SubElement(alignment, "align_2").text = result.add_gaps( gemmi.one_letter_code(chain.get_polymer().extract_sequence()), 2)
                iseq += 1


    return xml


def sequences(asuin):
    for seq in asuin.fileContent.seqList:
        residueKind, polymerType = {
            "PROTEIN": (gemmi.ResidueKind.AA, gemmi.PolymerType.PeptideL),
            "DNA": (gemmi.ResidueKind.DNA, gemmi.PolymerType.Dna),
            "RNA": (gemmi.ResidueKind.RNA, gemmi.PolymerType.Rna),
        }[seq.polymerType]
        expanded = gemmi.expand_one_letter_sequence(str(seq.sequence), residueKind)
        for _ in range(int(seq.nCopies)):
            yield expanded, polymerType
