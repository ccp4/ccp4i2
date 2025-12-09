"""
Split MTZ file into multiple mini-MTZ files based on column groups.

Uses gemmi directly via CCP4Utils.split_mtz_file for clean, CData-agnostic splitting.
"""
from __future__ import print_function

import os
import unittest
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4XtalData
from ccp4i2.core.CCP4Utils import split_mtz_file


# Map columnGroupType to output file class
COLUMN_GROUP_TYPE_MAP = {
    'Obs': CCP4XtalData.CObsDataFile,
    'Phs': CCP4XtalData.CPhsDataFile,
    'MapCoeffs': CCP4XtalData.CMapCoeffsDataFile,
    'FreeR': CCP4XtalData.CFreeRDataFile,
}

# Human-readable type descriptors for annotations
# Maps (group_type, subType) -> descriptor string
# subType values from CMapCoeffsDataFileStub: 1=normal, 2=difference, 3=anomalous difference
TYPE_DESCRIPTORS = {
    ('Obs', 1): 'Anom I',       # Anomalous intensities
    ('Obs', 2): 'Anom SF',      # Anomalous SFs
    ('Obs', 3): 'Mean I',       # Mean intensities
    ('Obs', 4): 'Mean SF',      # Mean SFs
    ('Phs', 1): 'HL Phs',       # Hendrickson-Lattman phases
    ('Phs', 2): 'Phi-FOM Phs',  # Phi/FOM phases
    ('MapCoeffs', 1): '2Fo-Fc',        # Normal map (e.g., FWT/PHWT)
    ('MapCoeffs', 2): 'Fo-Fc',         # Difference map (e.g., DELFWT/PHDELWT)
    ('MapCoeffs', 3): 'Anom diff',     # Anomalous difference map (e.g., FAN/PHAN)
    ('FreeR', None): 'FreeR',
}

# Column patterns for identifying map coefficient subtypes
# These are common naming conventions from refinement programs (refmac, phenix, buster)
MAP_SUBTYPE_PATTERNS = {
    # Difference map patterns (Fo-Fc / mFo-DFc) -> subType 2
    'difference': [
        'DELFWT', 'PHDELWT', 'DELPHWT',  # refmac standard
        'FOFC', 'FOFCWT', 'PHFOFCWT',     # phenix style
        'MFODFCWT', 'PHMFODFCWT',          # mFo-DFc weighted
        '2FOFCWT', 'PH2FOFCWT',            # 2Fo-Fc variants (actually normal maps)
    ],
    # Anomalous difference map patterns -> subType 3
    'anomalous': [
        'FAN', 'PHAN',           # anomalous difference
        'APTS', 'PHATS',         # anomalous phased target
        'APTS', 'SIGFAN',        # other anomalous
        'DANO', 'SIGDANO', 'PHIANO',  # anomalous difference
    ],
}


class splitMtz(CPluginScript):

    TASKTITLE = 'Import and Split MTZ to experimental data objects'
    TASKNAME = 'splitMtz'
    TASKCOMMAND = ''
    TASKVERSION = 0.1
    RUNEXTERNALPROCESS = False
    MAINTAINER = 'liz.potterton@york.ac.uk'

    ERROR_CODES = {
        200: {'description': 'No column groups selected for splitting'},
        201: {'description': 'Unknown column group type'},
        202: {'description': 'No columns found in column group'},
        203: {'description': 'Failed to split MTZ file'},
        204: {'description': 'Input MTZ file not found'},
    }

    def startProcess(self, command, **kw):  # noqa: ARG002
        inp = self.container.inputData
        out = self.container.outputData

        # Validate input file
        input_mtz = str(inp.HKLIN.fullPath)
        if not Path(input_mtz).exists():
            self.appendErrorReport(204, 'Input file: ' + input_mtz)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        input_basename = Path(str(inp.HKLIN.baseName)).stem
        files_created = 0

        # Process each selected column group
        for idx, group in enumerate(inp.COLUMNGROUPLIST):
            # Check if this group is selected
            selected = group.selected
            if hasattr(selected, 'value'):
                selected = selected.value
            if hasattr(selected, '__bool__'):
                selected = bool(selected)

            if not selected:
                continue

            # Get column group type and corresponding output class
            group_type = str(group.columnGroupType)
            output_cls = COLUMN_GROUP_TYPE_MAP.get(group_type)

            if output_cls is None:
                self.appendErrorReport(201, 'Unknown type: ' + group_type)
                continue

            # Extract column labels from the group
            column_labels = self._extract_column_labels(group)

            if not column_labels:
                self.appendErrorReport(202, 'Group index: ' + str(idx))
                continue

            # Split and create output file
            result = self._split_column_group(
                input_mtz, input_basename, group, column_labels, group_type, output_cls, out
            )
            if result == CPluginScript.FAILED:
                return CPluginScript.FAILED

            files_created += 1

        # Check we produced at least one output
        if files_created == 0:
            self.appendErrorReport(200, 'No column groups were marked as selected')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        print('Successfully created ' + str(files_created) + ' output file(s)')
        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def _extract_column_labels(self, group):
        """Extract column labels from a column group."""
        column_labels = []
        for col in group.columnList:
            label = col.columnLabel
            if hasattr(label, 'value'):
                label = label.value
            column_labels.append(str(label))
        return column_labels

    def _infer_map_subtype(self, column_labels):
        """
        Infer map coefficient subType from column labels.

        Returns:
            int: 1 = normal map (2Fo-Fc), 2 = difference map (Fo-Fc), 3 = anomalous difference
        """
        # Check column labels against known patterns (case-insensitive)
        upper_labels = [label.upper() for label in column_labels]

        # Check for anomalous patterns first (more specific)
        for label in upper_labels:
            for pattern in MAP_SUBTYPE_PATTERNS['anomalous']:
                if pattern.upper() in label:
                    return CCP4XtalData.CMapCoeffsDataFile.SUBTYPE_ANOM_DIFFERENCE  # 3

        # Check for difference map patterns
        for label in upper_labels:
            for pattern in MAP_SUBTYPE_PATTERNS['difference']:
                if pattern.upper() in label:
                    # Special case: 2FOFCWT is actually a normal map, not difference
                    if '2FOFC' in label or '2FO-FC' in label:
                        return CCP4XtalData.CMapCoeffsDataFile.SUBTYPE_NORMAL  # 1
                    return CCP4XtalData.CMapCoeffsDataFile.SUBTYPE_DIFFERENCE  # 2

        # Default to normal map (2Fo-Fc style)
        return CCP4XtalData.CMapCoeffsDataFile.SUBTYPE_NORMAL  # 1

    def _get_type_descriptor(self, group_type, content_flag, column_labels):
        """
        Get human-readable type descriptor for annotation.

        For MapCoeffs, infers subType from column names.
        For other types, uses contentFlag.
        """
        if group_type == 'MapCoeffs':
            subtype = self._infer_map_subtype(column_labels)
            return TYPE_DESCRIPTORS.get((group_type, subtype), group_type), subtype
        elif group_type == 'FreeR':
            return TYPE_DESCRIPTORS.get((group_type, None), group_type), None
        else:
            # For Obs and Phs, use contentFlag directly
            return TYPE_DESCRIPTORS.get((group_type, content_flag), group_type), None

    def _split_column_group(self, input_mtz, input_basename, group, column_labels, group_type, output_cls, out):
        """Split a single column group and create output file."""
        # Build output filename from dataset and column names
        dataset = str(group.dataset) if group.dataset else 'data'
        col_suffix = '_'.join(column_labels)
        output_filename = dataset + '_' + col_suffix + '.mtz'
        # Clean up any parentheses in filename
        output_filename = output_filename.replace('(', '').replace(')', '')
        output_path = os.path.join(self.workDirectory, output_filename)

        # Create column mapping (keep same names)
        column_mapping = {col: col for col in column_labels}

        # Get content flag
        content_flag = group.contentFlag
        if hasattr(content_flag, 'value'):
            content_flag = content_flag.value
        if hasattr(content_flag, '__int__'):
            content_flag = int(content_flag)

        # Get type descriptor and inferred subtype (for MapCoeffs)
        type_descriptor, inferred_subtype = self._get_type_descriptor(
            group_type, content_flag, column_labels
        )

        print('Splitting columns ' + str(column_labels) + ' to ' + output_path)

        try:
            # Use gemmi to split the MTZ
            split_mtz_file(input_mtz, output_path, column_mapping)

            # Create output file object and add to output list
            output_file = output_cls(fullPath=output_path)
            out.MINIMTZOUTLIST.append(output_file)

            # Set metadata on the output file
            out.MINIMTZOUTLIST[-1].contentFlag.set(content_flag)

            # Set subType for MapCoeffs based on column name heuristics
            if group_type == 'MapCoeffs' and inferred_subtype is not None:
                out.MINIMTZOUTLIST[-1].subType.set(inferred_subtype)

            # Annotation format: {type_descriptor} from {input_file} columns {dataset}/[col1,col2,...]
            col_list_str = ','.join(column_labels)
            annotation = type_descriptor + ' from ' + input_basename + ' columns ' + dataset + '/[' + col_list_str + ']'
            out.MINIMTZOUTLIST[-1].annotation.set(annotation)

            print('Created: ' + output_path + ' (' + type_descriptor + ')')
            return CPluginScript.SUCCEEDED

        except Exception as e:
            self.appendErrorReport(203, 'Error splitting ' + str(column_labels) + ': ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED


# ====================================================================================================
# PLUGIN TESTS


class testsplitMtz(unittest.TestCase):

    def setUp(self):
        from ccp4i2.core.CCP4Modules import QTAPPLICATION, PROCESSMANAGER
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)

    def tearDown(self):
        from ccp4i2.core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().setWaitForFinished(-1)


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testsplitMtz)
    return suite


def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
