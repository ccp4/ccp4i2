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
import xml.etree.ElementTree as ET

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

        xml = ET.Element("modelASUCheck")
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

    xml = ET.Element("SequenceAlignment")
    for irow, icol in zip(row_ind, col_ind):
        result, seq, chain = align_results[irow][icol]

        alignment = ET.SubElement(xml, "Alignment")
        ET.SubElement(alignment, "ChainID").text = chain.name
        ET.SubElement(alignment, "match_count").text = str(result.match_count)
        ET.SubElement(alignment, "identity").text = str(result.calculate_identity())
        ET.SubElement(alignment, "identity_1").text = str(result.calculate_identity(1))
        ET.SubElement(alignment, "identity_2").text = str(result.calculate_identity(2))
        ET.SubElement(alignment, "CIGAR").text = str(result.cigar_str())
        ET.SubElement(alignment, "align_1").text = result.add_gaps(
            gemmi.one_letter_code(seq), 1
        )
        ET.SubElement(alignment, "align_match").text = result.match_string
        ET.SubElement(alignment, "align_2").text = result.add_gaps(
            gemmi.one_letter_code(chain.get_polymer().extract_sequence()), 2
        )
    return xml


def sequences(asuin):
    for seq in asuin.fileContent.seqList:
        for _ in range(int(seq.nCopies)):
            if seq.polymerType == "PROTEIN":
                expanded = gemmi.expand_one_letter_sequence(seq.sequence, gemmi.ResidueKind.AA)
                yield expanded, gemmi.PolymerType.PeptideL
            elif seq.polymerType == "DNA":
                expanded = gemmi.expand_one_letter_sequence(seq.sequence, gemmi.ResidueKind.DNA)
                yield expanded, gemmi.PolymerType.Dna
            elif seq.polymerType == "RNA":
                expanded = gemmi.expand_one_letter_sequence(seq.sequence, gemmi.ResidueKind.RNA)
                yield expanded, gemmi.PolymerType.Rna
