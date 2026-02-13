"""
Implementation classes for CCP4ModelData.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from typing import Optional

from ccp4i2.core.cdata_stubs.CCP4ModelData import CAsuContentStub, CAsuContentSeqStub, CAsuContentSeqListStub, CAsuDataFileStub, CAtomRefmacSelectionStub, CAtomRefmacSelectionGroupsStub, CAtomRefmacSelectionListStub, CAtomRefmacSelectionOccupancyStub, CAtomSelectionStub, CBlastDataStub, CBlastDataFileStub, CBlastItemStub, CChemCompStub, CDictDataStub, CDictDataFileStub, CElementStub, CEnsembleStub, CEnsembleListStub, CEnsemblePdbDataFileStub, CHhpredDataStub, CHhpredDataFileStub, CHhpredItemStub, CMDLMolDataFileStub, CMol2DataFileStub, CMonomerStub, COccRefmacSelectionListStub, COccRelationRefmacListStub, CPdbDataStub, CPdbDataFileStub, CPdbDataFileListStub, CPdbEnsembleItemStub, CResidueRangeStub, CResidueRangeListStub, CSeqAlignDataFileStub, CSeqDataFileStub, CSeqDataFileListStub, CSequenceStub, CSequenceAlignmentStub, CSequenceMetaStub, CSequenceStringStub, CTLSDataFileStub


class CAsuContent(CAsuContentStub):
    """
    Base class for classes holding file contents

    Extends CAsuContentStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def loadFile(self, file_path: str = None):
        """
        Load ASU XML file and populate seqList.

        This method:
        1. Reads ASU XML file using ElementTree
        2. Extracts sequence information
        3. Populates seqList attribute with CAsuContentSeq objects

        Args:
            file_path: Optional path to ASU XML file. If None, gets path from parent CDataFile.

        Returns:
            CErrorReport with any errors encountered

        Example:
            # Load from parent file's path
            >>> asu_file = CAsuDataFile()
            >>> asu_file.setFullPath('/path/to/sequences.asu.xml')
            >>> asu_file.loadFile()  # Base CDataFile.loadFile() calls content.loadFile()
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        from pathlib import Path

        error = CErrorReport()

        # If no path provided, get from parent CDataFile
        if file_path is None:
            parent = self.get_parent()
            if parent is not None and hasattr(parent, 'getFullPath'):
                file_path = parent.getFullPath()

        # Validate file path
        if not file_path:
            return error

        path_obj = Path(file_path)
        if not path_obj.exists() or not path_obj.is_file():
            error.append(
                klass=self.__class__.__name__,
                code=101,
                details=f"ASU XML file does not exist or is not a file: '{file_path}'",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            import xml.etree.ElementTree as ET

            tree = ET.parse(file_path)
            root = tree.getroot()

            # Parse sequences from XML
            # Handle namespace
            ns = {'ccp4': 'http://www.ccp4.ac.uk/ccp4ns'}
            # Check if root IS ccp4i2_body (common case), otherwise search for it
            if root.tag == 'ccp4i2_body':
                body = root
            else:
                body = root.find('.//ccp4i2_body', ns) or root.find('.//ccp4i2_body')

            if body is not None:
                seqList_elem = body.find('seqList')
                if seqList_elem is not None:
                    # Clear existing sequences before loading to avoid duplicates
                    self.seqList.clear()

                    for seq_elem in seqList_elem.findall('CAsuContentSeq'):
                        seq_obj = CAsuContentSeq(parent=self.seqList, name=None)

                        # Parse sequence fields using smart assignment
                        if seq_elem.find('sequence') is not None:
                            seq_obj.sequence = seq_elem.find('sequence').text or ''
                        if seq_elem.find('nCopies') is not None:
                            seq_obj.nCopies = int(seq_elem.find('nCopies').text or '1')
                        if seq_elem.find('polymerType') is not None:
                            seq_obj.polymerType = seq_elem.find('polymerType').text or ''
                        if seq_elem.find('name') is not None:
                            seq_obj.name = seq_elem.find('name').text or ''
                        if seq_elem.find('description') is not None:
                            seq_obj.description = seq_elem.find('description').text or ''

                        # Add to seqList
                        self.seqList.append(seq_obj)

        except ET.ParseError as e:
            error.append(
                klass=self.__class__.__name__,
                code=102,
                details=f"Error parsing ASU XML file '{file_path}': {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=103,
                details=f"Error reading ASU XML file '{file_path}': {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )

        return error


class CAsuContentSeq(CAsuContentSeqStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    Extends CAsuContentSeqStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def formattedSequence(self):
        """
        Format the sequence for display with line wrapping and spacing.

        Formats sequence with:
        - 60 characters per line (6 blocks of 10)
        - Space every 10 characters for readability

        Returns:
            str: Formatted sequence string
        """
        GROUP_SIZE = 10
        GROUPS_PER_LINE = 6

        # Get sequence value - handle both CString and plain string
        sequence = self.sequence
        if hasattr(sequence, 'value'):
            seq = str(sequence.value)
        else:
            seq = str(sequence)

        # Clean whitespace from sequence
        seq = seq.replace(' ', '').replace('\n', '').replace('\t', '')

        if not seq:
            return ''

        # Split into groups of 10
        groups = [seq[i:i + GROUP_SIZE] for i in range(0, len(seq), GROUP_SIZE)]

        # Build lines with 6 groups each
        lines = []
        for i in range(0, len(groups), GROUPS_PER_LINE):
            line_groups = groups[i:i + GROUPS_PER_LINE]
            lines.append(' '.join(line_groups))

        return '\n'.join(lines)

    def numberOfResidues(self, countMulti=False):
        """
        Calculate the number of residues in the sequence.

        Args:
            countMulti: If True, multiply by nCopies

        Returns:
            int: Number of residues (sequence length), optionally multiplied by nCopies
        """
        # Get sequence value - handle both CString and plain string
        sequence = self.sequence
        if hasattr(sequence, 'value'):
            seq_str = str(sequence.value)
        else:
            seq_str = str(sequence)

        # Count non-whitespace characters
        num_residues = len(seq_str.replace(' ', '').replace('\n', '').replace('\t', ''))

        # Multiply by nCopies if requested
        if countMulti and hasattr(self, 'nCopies') and self.nCopies.isSet():
            num_residues *= int(self.nCopies.value)

        return num_residues

    def validity(self):
        """
        Validate this ASU content sequence entry.

        Checks:
        - nCopies is defined and >= 1
        - polymerType is defined and valid (PROTEIN, DNA, or RNA)
        - name is not null or zero length
        - sequence can be parsed by BioPython given the polymer type

        Returns:
            CErrorReport containing validation errors/warnings
        """
        import re
        from ccp4i2.core.base_object.error_reporting import CErrorReport, SEVERITY_ERROR

        report = CErrorReport()
        valid_polymer_types = {'PROTEIN', 'DNA', 'RNA'}

        # Get full object path for error reporting (matches frontend _objectPath)
        base_path = self.object_path() if hasattr(self, 'object_path') else ''

        # Check nCopies >= 1
        n_copies = getattr(self, 'nCopies', None)
        if n_copies is not None and hasattr(n_copies, 'value') and hasattr(n_copies, 'isSet'):
            n_copies_val = int(n_copies.value) if n_copies.isSet() else 0
        elif n_copies is not None:
            try:
                n_copies_val = int(n_copies)
            except (ValueError, TypeError):
                n_copies_val = 0
        else:
            n_copies_val = 0

        if n_copies_val < 1:
            # Use field path for cell-level highlighting
            field_path = f"{base_path}.nCopies" if base_path else 'nCopies'
            report.append(
                klass=self.__class__.__name__,
                code=101,
                details='nCopies must be >= 1',
                name=field_path,
                severity=SEVERITY_ERROR
            )

        # Check polymerType is valid
        polymer_type = getattr(self, 'polymerType', None)
        if polymer_type is not None and hasattr(polymer_type, 'value'):
            pt_str = str(polymer_type.value).strip().upper()
        elif polymer_type is not None:
            pt_str = str(polymer_type).strip().upper()
        else:
            pt_str = ''

        if not pt_str or pt_str not in valid_polymer_types:
            field_path = f"{base_path}.polymerType" if base_path else 'polymerType'
            report.append(
                klass=self.__class__.__name__,
                code=102,
                details=f'Invalid or missing polymer type (must be PROTEIN, DNA, or RNA, got "{pt_str}")',
                name=field_path,
                severity=SEVERITY_ERROR
            )

        # Check name is not null or empty
        name = getattr(self, 'name', None)
        if name is not None and hasattr(name, 'value'):
            name_str = str(name.value).strip()
        elif name is not None:
            name_str = str(name).strip()
        else:
            name_str = ''

        if not name_str:
            field_path = f"{base_path}.name" if base_path else 'name'
            report.append(
                klass=self.__class__.__name__,
                code=103,
                details='Name is required',
                name=field_path,
                severity=SEVERITY_ERROR
            )

        # Check sequence can be parsed by BioPython
        sequence = getattr(self, 'sequence', None)
        if sequence is not None and hasattr(sequence, 'value'):
            seq_str = str(sequence.value).replace(' ', '').replace('\n', '').replace('\t', '').upper()
        elif sequence is not None:
            seq_str = str(sequence).replace(' ', '').replace('\n', '').replace('\t', '').upper()
        else:
            seq_str = ''

        if len(seq_str) <= 1:
            field_path = f"{base_path}.sequence" if base_path else 'sequence'
            report.append(
                klass=self.__class__.__name__,
                code=104,
                details='Sequence must have more than 1 residue',
                name=field_path,
                severity=SEVERITY_ERROR
            )
        elif pt_str in valid_polymer_types:
            # Validate sequence characters based on polymer type
            field_path = f"{base_path}.sequence" if base_path else 'sequence'
            if pt_str == 'PROTEIN':
                # Standard amino acid one-letter codes
                valid_chars = set('GALMFWKQESPVICYHRNDT')
                invalid_chars = set(seq_str) - valid_chars
                if invalid_chars:
                    report.append(
                        klass=self.__class__.__name__,
                        code=105,
                        details=f'Sequence contains invalid amino acid characters: {", ".join(sorted(invalid_chars))}',
                        name=field_path,
                        severity=SEVERITY_ERROR
                    )
            else:
                # DNA or RNA - valid nucleotide codes
                valid_chars = set('CAUGT')
                invalid_chars = set(seq_str) - valid_chars
                if invalid_chars:
                    report.append(
                        klass=self.__class__.__name__,
                        code=106,
                        details=f'Sequence contains invalid nucleotide characters: {", ".join(sorted(invalid_chars))}',
                        name=field_path,
                        severity=SEVERITY_ERROR
                    )

        return report

    def molecularWeight(self, polymerType="PROTEIN"):
        """
        Calculate molecular weight of the sequence.

        Args:
            polymerType: Type of polymer ("PROTEIN", "RNA", or "DNA")

        Returns:
            float: Molecular weight in Daltons, multiplied by nCopies

        Note:
            Uses BioPython's SeqUtils.molecular_weight() function.
            Filters out non-standard residues before calculation.
        """
        import re
        import Bio.SeqUtils

        # Get sequence value - handle both CString and plain string
        sequence = self.sequence
        if hasattr(sequence, 'value'):
            seq_str = str(sequence.value)
        else:
            seq_str = str(sequence)

        # Get nCopies value - handle both CInt and plain int
        n_copies = self.nCopies
        if hasattr(n_copies, 'value'):
            n_copies_val = float(n_copies.value)
        else:
            n_copies_val = float(n_copies)

        # Normalize polymerType - handle CString, plain string, and empty values
        if hasattr(polymerType, 'value'):
            polymer_type_str = str(polymerType.value)
        elif hasattr(polymerType, '__str__'):
            polymer_type_str = str(polymerType)
        else:
            polymer_type_str = polymerType

        # Default to PROTEIN if empty
        if not polymer_type_str or polymer_type_str.strip() == '':
            polymer_type_str = "PROTEIN"

        # BioPython is not robust to bad sequences - filter to valid residues only
        if polymer_type_str == "PROTEIN":
            # Standard amino acid codes
            seq = re.sub('[^GALMFWKQESPVICYHRNDT]', '', seq_str)
            wt = Bio.SeqUtils.molecular_weight(seq, seq_type='protein')
        else:
            # Nucleic acid codes
            seq = re.sub('[^CAUGT]', '', seq_str)
            wt = Bio.SeqUtils.molecular_weight(seq, seq_type=polymer_type_str)

        return wt * n_copies_val


class CAsuContentSeqList(CAsuContentSeqListStub):
    """
    A list with all items of one CData sub-class

    Extends CAsuContentSeqListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """Initialize CAsuContentSeqList with subItem qualifier."""
        super().__init__(parent=parent, name=name, **kwargs)
        # Set the subItem qualifier to specify what type of items this list contains
        # subItem must be a dict with 'class' key pointing to the class
        self.set_qualifier('subItem', {'class': CAsuContentSeq, 'qualifiers': {}})
        # Require at least one entry for ASU content
        self.set_qualifier('listMinLength', 1)

    # Note: validity() is inherited from CList which handles:
    # - listMinLength/listMaxLength checks
    # - Aggregating child validity errors

    def validityAsDict(self):
        """
        Validate the ASU content list and return result as a dict.

        This is a convenience method for the frontend that returns a dict format
        suitable for displaying validation status messages.

        Returns:
            dict: Validation result with keys:
                - valid (bool): True if all checks pass
                - message (str): Description of validation status
                - errors (list): List of specific validation errors
        """
        report = self.validity()

        if report.maxSeverity() >= 4:  # SEVERITY_ERROR = 4
            errors = [err.get('details', str(err)) for err in report._errors]
            return {
                'valid': False,
                'message': f'{len(errors)} validation error(s)',
                'errors': errors
            }

        return {
            'valid': True,
            'message': f'{len(self)} valid sequence(s)',
            'errors': []
        }

    def molecularWeight(self):
        """
        Calculate total molecular weight of all sequences in the list.

        Sums the molecular weight of each CAsuContentSeq item, using its
        polymerType to determine the calculation method.

        Returns:
            float: Total molecular weight in Daltons
        """
        total_weight = 0.0
        for seq_obj in self:
            # Get polymerType - handle CString wrapper
            polymer_type = seq_obj.polymerType
            if hasattr(polymer_type, 'value'):
                polymer_type = str(polymer_type.value)
            else:
                polymer_type = str(polymer_type) if polymer_type else "PROTEIN"
            total_weight += seq_obj.molecularWeight(polymer_type)
        return total_weight


class CAsuDataFile(CAsuDataFileStub):
    """
    A reference to an XML file with CCP4i2 Header

    Extends CAsuDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    # No custom overrides needed - base CDataFile.loadFile() handles everything!
    # Just need CAsuContent.loadFile() implementation below

    def _clear_metadata(self):
        """Clear the selection dict when the file is cleared.

        This override ensures that when a user clears an ASU data file reference,
        any per-sequence selection metadata is also cleared. This prevents stale
        selection state from persisting when a different file is later selected.
        """
        if hasattr(self, 'selection') and self.selection is not None:
            self.selection.unSet()

    def isSelected(self, seqObj) -> bool:
        """
        Check if a sequence object is selected.

        The selection is stored in self.selection as a dict mapping sequence names
        to boolean values. If selection is not set, all sequences are considered
        selected by default.

        Args:
            seqObj: A CAsuContentSeq object from fileContent.seqList

        Returns:
            True if the sequence is selected (or selection not set), False otherwise
        """
        # Get the sequence name
        name = str(seqObj.name) if hasattr(seqObj, 'name') else str(seqObj)

        # If selection is not set, all sequences are selected by default
        if not hasattr(self, 'selection') or not self.selection.isSet():
            return True

        # Look up in selection dict, default to True if not found
        return self.selection.get(name, True)

    def saveFile(self, dbInfo=None):
        """
        Save fileContent to an XML file.

        Args:
            dbInfo: Optional database info dict with keys: projectName, projectId, jobId, jobNumber
        """
        # Get the file path to write to
        file_path = self.getFullPath()
        if not file_path:
            print(f"Warning: Cannot save CAsuDataFile - no file path set")
            return

        # Ensure parent directory exists
        from pathlib import Path
        output_path = Path(file_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Build XML structure
        import xml.etree.ElementTree as ET
        from xml.dom import minidom

        # Create root with CCP4i2 header structure
        root = ET.Element('ccp4i2_body')
        root.set('id', 'CAsuDataFile')

        # Add seqList element
        seq_list_elem = ET.SubElement(root, 'seqList')

        # Export each sequence from fileContent.seqList
        if hasattr(self, 'fileContent') and self.fileContent is not None:
            if hasattr(self.fileContent, 'seqList') and self.fileContent.seqList is not None:
                for seq_obj in self.fileContent.seqList:
                    seq_elem = ET.SubElement(seq_list_elem, 'CAsuContentSeq')

                    # Export sequence fields
                    if hasattr(seq_obj, 'sequence') and seq_obj.isSet('sequence'):
                        seq_text = ET.SubElement(seq_elem, 'sequence')
                        seq_text.text = str(seq_obj.sequence.value if hasattr(seq_obj.sequence, 'value') else seq_obj.sequence)

                    # Use allowDefault=True to include default values like nCopies=1 and polymerType=PROTEIN
                    if hasattr(seq_obj, 'nCopies') and seq_obj.isSet('nCopies', allowDefault=True):
                        copies = ET.SubElement(seq_elem, 'nCopies')
                        copies.text = str(seq_obj.nCopies.value if hasattr(seq_obj.nCopies, 'value') else seq_obj.nCopies)

                    if hasattr(seq_obj, 'polymerType') and seq_obj.isSet('polymerType', allowDefault=True):
                        polymer = ET.SubElement(seq_elem, 'polymerType')
                        polymer.text = str(seq_obj.polymerType.value if hasattr(seq_obj.polymerType, 'value') else seq_obj.polymerType)

                    if hasattr(seq_obj, 'name') and seq_obj.isSet('name'):
                        name_elem = ET.SubElement(seq_elem, 'name')
                        name_elem.text = str(seq_obj.name.value if hasattr(seq_obj.name, 'value') else seq_obj.name)

                    if hasattr(seq_obj, 'description') and seq_obj.isSet('description'):
                        desc = ET.SubElement(seq_elem, 'description')
                        desc.text = str(seq_obj.description.value if hasattr(seq_obj.description, 'value') else seq_obj.description)

                    # Export source element if it exists
                    if hasattr(seq_obj, 'source') and seq_obj.isSet('source'):
                        source_elem = ET.SubElement(seq_elem, 'source')
                        if hasattr(seq_obj.source, 'baseName') and seq_obj.source.isSet('baseName'):
                            base = ET.SubElement(source_elem, 'baseName')
                            base.text = str(seq_obj.source.baseName.value if hasattr(seq_obj.source.baseName, 'value') else seq_obj.source.baseName)
                        if hasattr(seq_obj.source, 'relPath') and seq_obj.source.isSet('relPath'):
                            rel = ET.SubElement(source_elem, 'relPath')
                            rel.text = str(seq_obj.source.relPath.value if hasattr(seq_obj.source.relPath, 'value') else seq_obj.source.relPath)

        # Pretty print XML
        xml_string = ET.tostring(root, encoding='unicode')
        dom = minidom.parseString(xml_string)
        pretty_xml = dom.toprettyxml(indent="  ")

        # Write to file
        with open(file_path, 'w') as f:
            f.write(pretty_xml)

        print(f"Debug: Saved ASU data to {file_path}")

    def writeFasta(
        self,
        fileName: str,
        indx: int = -1,
        format: str = 'fasta',
        writeMulti: bool = False,
        polymerTypes: list = None
    ):
        """
        Write sequences to a FASTA or PIR format file.

        Args:
            fileName: Output file path
            indx: Index of specific sequence to write (-1 for all)
            format: Output format ('fasta' or 'pir')
            writeMulti: Write multiple copies based on nCopies
            polymerTypes: List of polymer types to include (default: PROTEIN, RNA, DNA)
        """
        if polymerTypes is None:
            polymerTypes = ["PROTEIN", "RNA", "DNA"]

        # Load the file to populate fileContent
        self.loadFile()

        # Check if fileContent was populated
        if not hasattr(self, 'fileContent') or self.fileContent is None:
            print(f"Warning: No fileContent loaded from {self.getFullPath()}")
            return

        # Get selection mode from qualifiers if available
        selectionMode = self.get_qualifier('selectionMode', default=0)

        text = ''

        if indx < 0:
            # Write all sequences
            if not hasattr(self.fileContent, 'seqList'):
                return  # No sequences to write

            for seqObj in self.fileContent.seqList:
                # Filter by polymer type
                if str(seqObj.polymerType) not in polymerTypes:
                    continue

                # Determine number of copies to write
                if writeMulti:
                    nCopies = int(seqObj.nCopies)
                else:
                    nCopies = min(1, int(seqObj.nCopies))

                # Write each copy
                for nC in range(nCopies):
                    name = str(seqObj.name)

                    # Check selection if available
                    if selectionMode == 0 or (not hasattr(self, 'selection')) or \
                       (not self.selection.isSet()) or self.selection.get(name, True):

                        # Write FASTA/PIR header
                        text += '>' + name + '\n'
                        if format == 'pir':
                            text += '\n'

                        # Write sequence in 60-character lines
                        seq = str(seqObj.sequence)
                        while len(seq) > 0:
                            text += seq[0:60] + '\n'
                            seq = seq[60:]
        else:
            # Write single sequence by index
            if not hasattr(self.fileContent, 'seqList') or \
               indx >= len(self.fileContent.seqList):
                return

            seqObj = self.fileContent.seqList[indx]

            # Write FASTA/PIR header
            text += '>' + str(seqObj.name) + '\n'
            if format == 'pir':
                text += '\n'

            # Write sequence in 60-character lines
            seq = str(seqObj.sequence)
            while len(seq) > 0:
                text += seq[0:60] + '\n'
                seq = seq[60:]

        # Save to file
        with open(fileName, 'w') as f:
            f.write(text)

    def writeArpPir(self, fileName: str, writeMulti: bool = False, polymerTypes: list = None):
        """
        Write sequences to ARP/wARP PIR format file.

        Legacy API compatibility method for arp_warp_classic plugin.
        This is a wrapper around writeFasta() with format='pir'.

        Args:
            fileName: Output file path
            writeMulti: Write multiple copies based on nCopies
            polymerTypes: List of polymer types to include (default: PROTEIN, RNA, DNA)
        """
        return self.writeFasta(fileName, indx=-1, format='pir', writeMulti=writeMulti, polymerTypes=polymerTypes)


class CAtomRefmacSelection(CAtomRefmacSelectionStub):
    """
    A residue range selection for rigid body groups
    
    Extends CAtomRefmacSelectionStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAtomRefmacSelectionGroups(CAtomRefmacSelectionGroupsStub):
    """
    A group selection for occupancy groups
    
    Extends CAtomRefmacSelectionGroupsStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAtomRefmacSelectionList(CAtomRefmacSelectionListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CAtomRefmacSelectionListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAtomRefmacSelectionOccupancy(CAtomRefmacSelectionOccupancyStub):
    """
    A residue range selection for occupancy groups
    
    Extends CAtomRefmacSelectionOccupancyStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAtomSelection(CAtomSelectionStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CAtomSelectionStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CBlastData(CBlastDataStub):
    """
    Base class for classes holding file contents
    
    Extends CBlastDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CBlastDataFile(CBlastDataFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CBlastDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CBlastItem(CBlastItemStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CBlastItemStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CChemComp(CChemCompStub):
    """
    Component of CDictDataFile contents
    
    Extends CChemCompStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDictData(CDictDataStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CDictDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDictDataFile(CDictDataFileStub):
    """
    A refmac dictionary file
    
    Extends CDictDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CElement(CElementStub):
    """
    Chemical element 
    
    Extends CElementStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CEnsemble(CEnsembleStub):
    """
    An ensemble of models. Typically, this would be a set of related
PDB files, but models could also be xtal or EM maps. This should
be indicated by the types entry.
A single ensemble is a CList of structures.

    Extends CEnsembleStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, *args, **kwargs):
        """Initialize CEnsemble with a default pdbItemList containing one item."""
        super().__init__(*args, **kwargs)

        # Set the subItem qualifier for pdbItemList to specify that it contains CPdbEnsembleItem objects
        if self.pdbItemList is not None:
            self.pdbItemList.set_qualifier('subItem', {
                'class': CPdbEnsembleItem,  # Use implementation class, not stub
                'qualifiers': {}
            })

        # Set default number of copies to 1 (not 0, which means "no search")
        if hasattr(self, 'number') and self.number is not None:
            # Check if number is unset OR explicitly set to 0
            if not self.number.isSet() or self.number.value == 0:
                self.number.value = 1

            # NOTE: Temporarily disabled default item creation to fix i2run nested path parsing
            # Legacy code expects pdbItemList to have at least one item by default
            # phaser_simple.py line 33-35 assumes: elements = ensemble.pdbItemList; pdbItem = elements[-1]
            # However, this causes issues with i2run argument parsing creating duplicate items
            # The default item gets created, i2run sets attributes on it, but then XML shows both
            # the original empty state AND the modified state as separate items
            # TODO: Investigate why XML serialization duplicates the item
            # if len(self.pdbItemList) == 0:
            #     # Add one default item to the pdbItemList
            #     self.pdbItemList.append(self.pdbItemList.makeItem())


class CEnsembleList(CEnsembleListStub):
    """
    A list with all items of one CData sub-class

    Extends CEnsembleListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, *args, **kwargs):
        """Initialize CEnsembleList with CEnsemble as the subItem type."""
        super().__init__(*args, **kwargs)

        # Set the subItem qualifier to specify that this list contains CEnsemble objects
        # This is required for makeItem() to create the correct type
        self.set_qualifier('subItem', {
            'class': CEnsemble,
            'qualifiers': {}
        })

    def makeItem(self, *args, **kwargs):
        """
        Create a new CEnsemble item with auto-generated label.

        Legacy expectation: When a CEnsemble is created via makeItem(),
        it should automatically get a label "Ensemble_{index}" where
        index is the current length of the list.

        Returns:
            CEnsemble: A new ensemble with auto-generated label
        """
        # Create the item using parent's makeItem
        item = super().makeItem(*args, **kwargs)

        # Auto-generate label based on current list length
        if hasattr(item, 'label') and item.label is not None:
            if not item.label.isSet():
                index = len(self)  # Current length before appending
                item.label.value = f"Ensemble_{index}"

        return item


class CEnsemblePdbDataFile(CEnsemblePdbDataFileStub):
    """
    A PDB coordinate file containing ensemble of structures as 'NMR' models
    
    Extends CEnsemblePdbDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CHhpredData(CHhpredDataStub):
    """
    Base class for classes holding file contents
    
    Extends CHhpredDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CHhpredDataFile(CHhpredDataFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CHhpredDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CHhpredItem(CHhpredItemStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CHhpredItemStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMDLMolDataFile(CMDLMolDataFileStub):
    """
    A molecule definition file (MDL)

    Extends CMDLMolDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set default file extension for MDL MOL files
        self.set_qualifier('fileExtensions', ['mol'])
        self.set_qualifier('mimeTypeName', 'chemical/x-mdl-molfile')


class CMol2DataFile(CMol2DataFileStub):
    """
    A molecule definition file (MOL2)

    Extends CMol2DataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set default file extension for MOL2 files
        self.set_qualifier('fileExtensions', ['mol2'])
        self.set_qualifier('mimeTypeName', 'chemical/x-mol2')


class CMonomer(CMonomerStub):
    """
    A monomer compound. ?smiles
    
    Extends CMonomerStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class COccRefmacSelectionList(COccRefmacSelectionListStub):
    """
    A list with all items of one CData sub-class
    
    Extends COccRefmacSelectionListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class COccRelationRefmacList(COccRelationRefmacListStub):
    """
    A list with all items of one CData sub-class
    
    Extends COccRelationRefmacListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


# Try to import BioPython - it's optional but provides enhanced functionality
try:
    import Bio.SeqIO
    import Bio.AlignIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

# BioPython format mappings
SEQFORMATLIST = ['fasta', 'pir', 'swiss', 'tab', 'ig']
ALIGNFORMATLIST = ['fasta', 'clustal', 'phylip', 'stockholm', 'nexus', 'emboss', 'msf']
EXTLIST = {
    'fasta': 'fasta', 'fa': 'fasta', 'faa': 'fasta', 'fas': 'fasta',
    'pir': 'pir',
    'aln': 'clustal', 'clw': 'clustal', 'clustal': 'clustal',
    'phy': 'phylip', 'phylip': 'phylip',
    'stockholm': 'stockholm', 'stk': 'stockholm',
    'nexus': 'nexus', 'nex': 'nexus',
    'msf': 'msf'
}


class CBioPythonSeqInterface:
    """
    Mixin class providing BioPython sequence loading functionality.

    This interface allows loading sequences from various file formats using BioPython.
    It provides format detection, file fixing for malformed files, and robust error handling.
    """

    def simpleFormatTest(self, filename: str) -> str:
        """
        Detect sequence file format by inspecting content.

        Args:
            filename: Path to sequence file

        Returns:
            Format string: 'fasta', 'pir', 'noformat', or 'unknown'
        """
        import re
        from pathlib import Path

        try:
            text = Path(filename).read_text()
        except Exception:
            return 'unknown'

        segments = text.split('>')
        if len(segments[0]) != 0:
            return 'unknown'

        lines = ('>' + segments[1]).split('\n')
        nonAlphaList = []

        for il, line in enumerate(lines):
            # If only non-alpha character is a "-" and it's not the first line, treat as gap
            if il > 0 and ((len(set(re.findall('[^(a-z,A-Z, )]', line))) == 1 and
                           list(set(re.findall('[^(a-z,A-Z, )]', line)))[0] == "-") or
                          line.endswith("*")):
                nonAlphaList.append(0)
            else:
                nonAlphaList.append(len(re.findall('[^(a-z,A-Z, )]', line)))

        nonAlphaTot = sum(nonAlphaList)

        if nonAlphaTot == 0:
            return 'noformat'
        elif len(lines) > 0 and len(lines[0]) > 0 and lines[0][0] == '>':
            if ';' in lines[0] and (nonAlphaTot - (nonAlphaList[0] + nonAlphaList[1] + nonAlphaList[-1])) == 0:
                return 'pir'
            elif (nonAlphaTot - nonAlphaList[0]) == 0:
                return 'fasta'

        return 'unknown'

    def fixPirFile(self, filename: str):
        """
        Fix malformed PIR format files.

        PIR files should have format:
        >P1;identifier
        description
        sequence...
        *

        Args:
            filename: Path to PIR file to fix

        Returns:
            Tuple of (fixed_file_path, CErrorReport, data_dict) or (None, error, None)
        """
        import tempfile
        import os
        from pathlib import Path
        from ccp4i2.core.base_object.error_reporting import CErrorReport

        try:
            text = Path(filename).read_text()
        except Exception:
            return None, None, None

        if len(text.strip()) == 0:
            return None, CErrorReport(klass='CBioPythonSeqInterface', code=412,
                                     details='Sequence file is empty'), None

        fragments = text.split('\n>')
        if len(fragments) > 0 and len(fragments[0]) > 0 and fragments[0][0] == '>':
            fragments[0] = fragments[0][1:]

        output = ''
        for text_frag in fragments:
            text_frag = text_frag.strip()
            if len(text_frag) > 2 and text_frag[2] != ';':
                text_frag = 'P1;' + text_frag
            if '*' not in text_frag:
                text_frag = text_frag + '*'
            output = output + '>' + text_frag + '\n'

        # Write to temporary file
        fd, temp_path = tempfile.mkstemp(suffix='.pir')
        try:
            os.write(fd, output.encode('utf-8'))
            os.close(fd)

            err, data = self.bioLoadSeqFile(temp_path, 'pir')
            if err.count() == 0:
                return temp_path, err, data
            else:
                os.unlink(temp_path)
                return None, None, None
        except Exception:
            if fd:
                os.close(fd)
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            return None, None, None

    def fixFastaFile(self, filename: str):
        """
        Fix plain sequence files by adding FASTA header.

        Args:
            filename: Path to plain sequence file

        Returns:
            Tuple of (fixed_file_path, CErrorReport, data_dict) or (None, error, None)
        """
        import tempfile
        import os
        import re
        from pathlib import Path
        from ccp4i2.core.base_object.error_reporting import CErrorReport

        try:
            text = Path(filename).read_text()
        except Exception:
            return None, None, None

        if len(text.strip()) == 0:
            return None, CErrorReport(klass='CBioPythonSeqInterface', code=412,
                                     details='Sequence file is empty'), None

        lines = text.split('\n')
        nonAlphaList = []
        for line in lines:
            nonAlphaList.extend(re.findall('[^(a-z,A-Z, )]', line))

        if len(nonAlphaList) > 0:
            return None, None, None

        # Add FASTA header using filename as identifier
        output = '>' + Path(filename).stem + '\n' + text

        # Write to temporary file
        fd, temp_path = tempfile.mkstemp(suffix='.fasta')
        try:
            os.write(fd, output.encode('utf-8'))
            os.close(fd)

            err, data = self.bioLoadSeqFile(temp_path, 'fasta')
            if err.count() == 0:
                return temp_path, err, data
            else:
                os.unlink(temp_path)
                return None, None, None
        except Exception:
            if fd:
                os.close(fd)
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            return None, None, None

    def bioLoadSeqFile(self, filename: str, format: str, record: int = 0):
        """
        Load sequence file using BioPython.

        Args:
            filename: Path to sequence file
            format: BioPython format string ('fasta', 'pir', 'clustal', etc.)
            record: Index of record to load (for multi-record files)

        Returns:
            Tuple of (CErrorReport, data_dict) where data_dict contains:
                - format: Format used
                - identifier: Sequence identifier
                - sequence: Sequence string
                - name: Sequence name
                - description: Sequence description
                - idList: List of all identifiers in file
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport

        err = CErrorReport()

        if not BIOPYTHON_AVAILABLE:
            err.append(klass='CBioPythonSeqInterface', code=402,
                      details='BioPython not available - cannot load sequence file')
            return err, {}

        try:
            with open(filename, 'r') as f:
                if format in SEQFORMATLIST:
                    seq_records = list(Bio.SeqIO.parse(f, format))
                elif format in ALIGNFORMATLIST:
                    ali_records = list(Bio.AlignIO.parse(f, format))
                    if len(ali_records) > 0:
                        seq_records = list(ali_records[0])
                    else:
                        seq_records = []
                else:
                    err.append(klass='CBioPythonSeqInterface', code=403,
                              details=f'Unknown sequence file format: {format}')
                    return err, {}
        except Exception as e:
            err.append(klass='CBioPythonSeqInterface', code=402,
                      details=f'Error reading sequence file: {e}')
            return err, {}

        if len(seq_records) == 0:
            err.append(klass='CBioPythonSeqInterface', code=412,
                      details='Sequence file is empty or contains no valid records')
            return err, {}

        if record >= len(seq_records):
            err.append(klass='CBioPythonSeqInterface', code=405,
                      details=f'Record index {record} out of range (file has {len(seq_records)} records)')
            return err, {}

        idList = [rec.id for rec in seq_records]
        selected_record = seq_records[record]

        return err, {
            'format': format,
            'identifier': selected_record.id,
            'sequence': str(selected_record.seq),
            'name': selected_record.name,
            'description': selected_record.description,
            'idList': idList
        }

    def loadExternalFile(self, filename: str, format: str = None, record: int = 0):
        """
        Load sequence from external file with automatic format detection.

        Args:
            filename: Path to sequence file
            format: Format hint ('fasta', 'pir', 'unknown', etc.)
            record: Index of record to load (for multi-record files)

        Returns:
            Dictionary with sequence data or None on error
        """
        from pathlib import Path
        from ccp4i2.core.base_object.error_reporting import CErrorReport

        if not Path(filename).exists():
            return None

        # Guess format from extension if not provided
        if format is None or format == 'unknown':
            ext = Path(filename).suffix[1:].lower()
            format = EXTLIST.get(ext, None)

        # Build test order - try expected format first
        testOrder = SEQFORMATLIST + ALIGNFORMATLIST
        if format is not None and format in testOrder:
            testOrder.remove(format)
            testOrder.insert(0, format)

        # Try each format
        for test_format in testOrder:
            err, data = self.bioLoadSeqFile(filename, test_format, record=record)
            if err.count() == 0:
                return data

            # Special handling for PIR format
            if test_format == format and format == 'pir':
                fixed_path, fix_err, fix_data = self.fixPirFile(filename)
                if fixed_path is not None:
                    return fix_data

        # Try fixing as plain FASTA
        fixed_path, fix_err, fix_data = self.fixFastaFile(filename)
        if fixed_path is not None:
            return fix_data

        return None


class CPdbDataComposition:
    """
    Coordinate file composition analysis using gemmi.

    Analyzes PDB/mmCIF structure to extract:
    - chains: List of chain IDs
    - monomers: List of ligand residue IDs (chain:resname:seqnum format) for visualization
    - allResidueNames: List of all unique residue names (for backward compatibility)
    - peptides: List of chains containing amino acids
    - nucleics: List of chains containing nucleic acids
    - solventChains: List of chains containing solvent
    - saccharides: List of chains containing sugar residues
    - chainInfo: List of [nres, first_resid, last_resid] for each chain
    - nChains: Number of chains
    - nResidues: Total number of residues
    - nAtoms: Total number of atoms
    - nresSolvent: Number of solvent residues
    - elements: List of unique element types (excluding C, N, O)
    - containsHydrogen: Boolean for hydrogen presence

    The 'monomers' list now contains only true ligands:
    - Non-polymer entities (using gemmi's entity_type)
    - Not water molecules
    - Not single metal ions
    - With more than 5 atoms (significant ligands)

    Format: "chain:resname:seqnum" e.g. "A:ATP:501"
    """

    # Minimum atoms for a residue to be considered a "significant" ligand
    MIN_LIGAND_ATOMS = 5

    def __init__(self, gemmi_structure):
        """
        Analyze structure composition using gemmi.

        Args:
            gemmi_structure: gemmi.Structure object
        """
        import gemmi
        from ccp4i2.core.coordinate_selection.evaluator import (
            AMINO_ACIDS, NUCLEIC_ACIDS, SOLVENTS, SACCHARIDES
        )

        self.chains = []
        self.monomers = []  # Now: significant ligands only (chain:resname:seqnum format)
        self.allResidueNames = []  # All unique residue names (backward compat)
        self.peptides = []
        self.nucleics = []
        self.solventChains = []
        self.saccharides = []
        self.chainInfo = []
        self.nModels = len(gemmi_structure)
        self.nChains = 0
        self.nResidues = 0
        self.nAtoms = 0
        self.nresSolvent = 0
        self.elements = []
        self.containsHydrogen = False
        # Enhanced digest fields
        self.ligands = []           # [{chain, name, seqNum, atomCount}, ...]
        self.chainDetails = []      # [{id, type, nResidues, nAtoms, firstRes, lastRes, ligandCount, hasAltConf}, ...]
        self.residueNameCounts = {} # {resName: count, ...}

        all_resname_set = set()
        ligand_list = []
        structured_ligand_list = []
        residue_name_counts = {}
        element_set = set()

        # Analyze first model only (like legacy code)
        if len(gemmi_structure) > 0:
            model = gemmi_structure[0]
            self.nChains = len(model)

            for chain in model:
                chain_id = chain.name
                self.chains.append(chain_id)

                has_amino = False
                has_nucleic = False
                has_solvent = False
                has_saccharide = False
                has_altconf_chain = False
                ligand_count_chain = 0

                nres = 0
                first_resid = None
                last_resid = None
                nres_solvent_chain = 0
                natoms_chain = 0

                for residue in chain:
                    res_name = residue.name
                    all_resname_set.add(res_name)
                    residue_name_counts[res_name] = residue_name_counts.get(res_name, 0) + 1
                    nres += 1
                    natoms_res = len(residue)

                    # Track first and last residue IDs
                    resid_str = str(residue.seqid.num)
                    if first_resid is None:
                        first_resid = resid_str
                    last_resid = resid_str

                    # Use gemmi's entity_type for proper classification
                    entity_type = residue.entity_type

                    if entity_type == gemmi.EntityType.Polymer:
                        # Classify by residue name for accurate protein vs nucleic
                        if res_name in NUCLEIC_ACIDS:
                            has_nucleic = True
                        elif res_name in AMINO_ACIDS:
                            has_amino = True
                        else:
                            has_amino = True  # Fallback for unknown polymer
                    elif entity_type == gemmi.EntityType.Water:
                        has_solvent = True
                        nres_solvent_chain += 1
                    elif entity_type == gemmi.EntityType.NonPolymer:
                        # This is a ligand/hetero compound
                        # Check if it's a single metal ion
                        is_metal_ion = (natoms_res == 1 and
                                       len(residue) > 0 and
                                       residue[0].element.is_metal)

                        if not is_metal_ion and natoms_res >= self.MIN_LIGAND_ATOMS:
                            # Significant ligand - add to monomers list
                            ligand_id = f"{chain_id}:{res_name}:{residue.seqid.num}"
                            ligand_list.append(ligand_id)
                            structured_ligand_list.append({
                                "chain": chain_id,
                                "name": res_name,
                                "seqNum": residue.seqid.num,
                                "atomCount": natoms_res,
                            })
                            ligand_count_chain += 1

                        # Check for saccharides
                        if res_name in SACCHARIDES:
                            has_saccharide = True

                    # Count atoms and analyze elements
                    for atom in residue:
                        natoms_chain += 1
                        element = atom.element.name

                        # Track non-common elements
                        if element not in ['C', 'N', 'O']:
                            element_set.add(element)

                        # Check for hydrogen
                        if element in ['H', 'D']:
                            self.containsHydrogen = True

                        # Check for alternate conformations
                        if not has_altconf_chain and atom.altloc and atom.altloc not in ('\x00', ' ', ''):
                            has_altconf_chain = True

                # Store chain info
                self.chainInfo.append([nres, first_resid or '', last_resid or ''])
                self.nResidues += nres
                self.nAtoms += natoms_chain
                self.nresSolvent += nres_solvent_chain

                # Classify chains
                if has_amino:
                    self.peptides.append(chain_id)
                if has_nucleic:
                    self.nucleics.append(chain_id)
                if has_solvent:
                    self.solventChains.append(chain_id)
                if has_saccharide:
                    self.saccharides.append(chain_id)

                # Build per-chain detail
                chain_type = ("protein" if has_amino else
                              "nucleic" if has_nucleic else
                              "solvent" if has_solvent else
                              "saccharide" if has_saccharide else "other")
                self.chainDetails.append({
                    "id": chain_id,
                    "type": chain_type,
                    "nResidues": nres,
                    "nAtoms": natoms_chain,
                    "firstRes": first_resid or '',
                    "lastRes": last_resid or '',
                    "ligandCount": ligand_count_chain,
                    "hasAltConf": has_altconf_chain,
                })

        # monomers now contains only significant ligands in chain:resname:seqnum format
        self.monomers = sorted(ligand_list)
        self.allResidueNames = sorted(list(all_resname_set))
        self.elements = sorted(list(element_set))
        self.ligands = structured_ligand_list
        self.residueNameCounts = residue_name_counts


class CPdbData(CPdbDataStub):
    """
    Contents of a PDB file - a subset with functionality for GUI

    Extends CPdbDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def loadFile(self, file_path: str = None):
        """
        Load PDB or mmCIF coordinate file using gemmi library.

        This method:
        1. Reads coordinate file using gemmi.read_structure()
        2. Stores gemmi Structure object for queries
        3. Handles both PDB and mmCIF formats automatically

        Args:
            file_path: Optional path to coordinate file (.pdb, .cif, .ent). If None, gets path from parent CDataFile.

        Returns:
            CErrorReport with any errors encountered

        Example:
            # Load from parent file's path
            >>> pdb_file = CPdbDataFile()
            >>> pdb_file.setFullPath('/path/to/model.pdb')
            >>> pdb_file.fileContent.loadFile()

            # Load from explicit path (legacy pattern)
            >>> pdb_data = CPdbData()
            >>> error = pdb_data.loadFile('/path/to/model.pdb')
            >>> if error.count() == 0:
            ...     print(f"Loaded structure with {len(pdb_data._gemmi_structure)} models")
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        from pathlib import Path

        error = CErrorReport()

        # If no path provided, get from parent CDataFile
        if file_path is None:
            parent = self.get_parent()
            if parent is not None and hasattr(parent, 'getFullPath'):
                file_path = parent.getFullPath()

        # Validate file path
        if not file_path:
            return error

        path_obj = Path(file_path)
        if not path_obj.exists() or not path_obj.is_file():
            error.append(
                klass=self.__class__.__name__,
                code=101,
                details=f"Coordinate file does not exist or is not a file: '{file_path}'",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            import gemmi
        except ImportError as e:
            error.append(
                klass=self.__class__.__name__,
                code=102,
                details=f"Failed to import gemmi library: {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            # Read structure using gemmi (handles PDB and mmCIF automatically)
            structure = gemmi.read_structure(str(file_path))

            # Store gemmi Structure object for advanced queries
            # Use object.__setattr__ to bypass smart assignment
            object.__setattr__(self, '_gemmi_structure', structure)

            # Create composition analysis
            object.__setattr__(self, '_composition', CPdbDataComposition(structure))

            # Emit signal if available
            if hasattr(self, 'dataChanged'):
                try:
                    self.dataChanged.emit()
                except:
                    pass  # Signal system may not be available

        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=103,
                details=f"Error reading coordinate file '{file_path}': {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )

        return error

    @property
    def composition(self):
        """
        Get composition analysis of the coordinate file.

        Returns:
            CPdbDataComposition instance with chains, monomers, etc., or None if not loaded
        """
        return getattr(self, '_composition', None)

    def isMMCIF(self) -> bool:
        """
        Check if the loaded coordinate file is in mmCIF format.

        This method checks the gemmi structure's input format to determine
        if the file was mmCIF (vs PDB).

        Returns:
            bool: True if file is mmCIF format, False otherwise (including if not loaded)

        Example:
            >>> pdb_data = CPdbData()
            >>> pdb_data.loadFile('/path/to/model.cif')
            >>> if pdb_data.isMMCIF():
            ...     print("File is in mmCIF format")
        """
        # Check if gemmi structure is loaded
        gemmi_structure = getattr(self, '_gemmi_structure', None)
        if gemmi_structure is None:
            return False

        # Try to determine format from gemmi structure metadata
        # gemmi.Structure has an 'input_format' attribute that indicates the original format
        if hasattr(gemmi_structure, 'input_format'):
            import gemmi
            # gemmi.CoorFormat.Mmcif indicates mmCIF format
            return gemmi_structure.input_format == gemmi.CoorFormat.Mmcif

        # Fallback: check file name if parent file exists
        parent = self.get_parent()
        if parent is not None and hasattr(parent, 'getFullPath'):
            file_path = parent.getFullPath()
            if file_path:
                from pathlib import Path
                suffix = Path(file_path).suffix.lower()
                return suffix in ['.cif', '.mmcif']

        return False

    def interpretSelection(self, selection_string: str):
        """
        Parse and apply a coordinate selection string.

        This method parses an mmdb-style CID selection string and returns
        the selected atoms from the loaded structure.

        Args:
            selection_string: Selection string (e.g., "A/27.A", "{A/ or B/} and {(ALA)}")

        Returns:
            Tuple of (num_atoms, selected_atoms_list)
            where selected_atoms_list is list of (model, chain, residue, atom) tuples

        Raises:
            ValueError: If selection string is invalid or structure not loaded

        Example:
            >>> pdb_data = CPdbData()
            >>> pdb_data.loadFile('/path/to/model.pdb')
            >>> n_atoms, atoms = pdb_data.interpretSelection("A/27.A")
            >>> print(f"Selected {n_atoms} atoms")
        """
        # Ensure structure is loaded
        gemmi_structure = getattr(self, '_gemmi_structure', None)
        if gemmi_structure is None:
            raise ValueError("Structure not loaded - call loadFile() first")

        try:
            # Parse and evaluate selection using our coordinate_selection module
            from ccp4i2.core.coordinate_selection import parse_selection, evaluate_selection

            ast = parse_selection(selection_string)
            selected_atoms = evaluate_selection(ast, gemmi_structure)

            return (len(selected_atoms), selected_atoms)

        except Exception as e:
            raise ValueError(f"Failed to interpret selection '{selection_string}': {e}")

    def writeSelection(self, selected_atoms, file_path: str):
        """
        Write selected atoms to a PDB or mmCIF file.

        Args:
            selected_atoms: List of (model, chain, residue, atom) tuples from interpretSelection()
            file_path: Output file path (.pdb or .cif extension)

        Returns:
            0 on success, non-zero on failure

        Example:
            >>> n_atoms, atoms = pdb_data.interpretSelection("A/")
            >>> pdb_data.writeSelection(atoms, "/tmp/chain_A.pdb")
        """
        import gemmi
        from pathlib import Path

        try:
            # Create a new structure for output
            output_structure = gemmi.Structure()

            # Copy metadata from original structure
            gemmi_structure = getattr(self, '_gemmi_structure', None)
            if gemmi_structure:
                output_structure.name = gemmi_structure.name
                output_structure.cell = gemmi_structure.cell
                output_structure.spacegroup_hm = gemmi_structure.spacegroup_hm

            # Group atoms by model and chain for efficient writing
            from collections import defaultdict
            by_model_chain = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

            for model, chain, residue, atom in selected_atoms:
                model_idx = 0  # Use first model for simplicity
                for idx, m in enumerate(gemmi_structure):
                    if m == model:
                        model_idx = idx
                        break

                chain_id = chain.name
                res_key = (residue.seqid.num, residue.seqid.icode, residue.name)
                by_model_chain[model_idx][chain_id][res_key].append(atom)

            # Build output structure
            for model_idx in sorted(by_model_chain.keys()):
                output_model = output_structure.add_model(gemmi.Model(str(model_idx + 1)))

                for chain_id in sorted(by_model_chain[model_idx].keys()):
                    output_chain = output_model.add_chain(gemmi.Chain(chain_id))

                    for res_key in sorted(by_model_chain[model_idx][chain_id].keys()):
                        seq_num, ins_code, res_name = res_key
                        seqid = gemmi.SeqId(seq_num, ins_code if ins_code else ' ')
                        output_residue = output_chain.add_residue(gemmi.Residue())
                        output_residue.name = res_name
                        output_residue.seqid = seqid

                        for atom in by_model_chain[model_idx][chain_id][res_key]:
                            # Create a copy of the atom
                            new_atom = gemmi.Atom()
                            new_atom.name = atom.name
                            new_atom.element = atom.element
                            new_atom.pos = atom.pos
                            new_atom.occ = atom.occ
                            new_atom.b_iso = atom.b_iso
                            new_atom.charge = atom.charge
                            new_atom.altloc = atom.altloc
                            output_residue.add_atom(new_atom)

            # Determine output format from extension
            suffix = Path(file_path).suffix.lower()
            if suffix in ['.cif', '.mmcif']:
                output_structure.make_mmcif_document().write_file(str(file_path))
            else:
                # Default to PDB format
                output_structure.write_pdb(str(file_path))

            return 0

        except Exception as e:
            print(f"Error writing selection to {file_path}: {e}")
            import traceback
            traceback.print_exc()
            return 1


class CPdbDataFile(CPdbDataFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    Extends CPdbDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: fileContentClassName='CPdbData' is already set in CPdbDataFileStub decorator

    def isSelectionSet(self) -> bool:
        """Check if an atom selection is defined for this PDB file.

        Returns True if the selection.text attribute is set and contains
        non-whitespace content after stripping leading/trailing whitespace.

        Returns:
            bool: True if selection is set and non-empty, False otherwise
        """
        import sys
        # Check if selection attribute exists
        if not hasattr(self, 'selection'):
            pass  # DEBUG: print(f"DEBUG isSelectionSet: No selection attribute", file=sys.stderr)
            return False

        # Check if text attribute exists
        if not hasattr(self.selection, 'text'):
            pass  # DEBUG: print(f"DEBUG isSelectionSet: selection has no text attribute", file=sys.stderr)
            return False

        # Get the text value
        text_value = self.selection.text.value
        if text_value is None:
            pass  # DEBUG: print(f"DEBUG isSelectionSet: text_value is None", file=sys.stderr)
            return False

        # Check if text has non-whitespace content
        stripped = str(text_value).strip()
        result = len(stripped) > 0
        pass  # DEBUG: print(f"DEBUG isSelectionSet: text_value='{text_value}', stripped='{stripped}', result={result}", file=sys.stderr)
        return result

    def _introspect_content_flag(self) -> Optional[int]:
        """Auto-detect contentFlag by determining if file is PDB or mmCIF format.

        Returns:
            1 (CONTENT_FLAG_PDB) if PDB format
            2 (CONTENT_FLAG_MMCIF) if mmCIF format
            None if file cannot be read or format cannot be determined
        """
        from pathlib import Path

        file_path = self.getFullPath()
        if not file_path or not Path(file_path).exists():
            return None

        try:
            # Check the file extension first
            path = Path(file_path)
            suffix = path.suffix.lower()

            # mmCIF files typically have .cif or .mmcif extensions
            if suffix in ['.cif', '.mmcif']:
                return self.__class__.CONTENT_FLAG_MMCIF

            # PDB files typically have .pdb or .ent extensions
            elif suffix in ['.pdb', '.ent']:
                return self.__class__.CONTENT_FLAG_PDB

            # If extension is ambiguous, check file content
            # mmCIF files start with "data_" or have loop structures
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('data_') or first_line.startswith('loop_'):
                    return self.__class__.CONTENT_FLAG_MMCIF
                # PDB files typically start with record types like HEADER, CRYST1, ATOM, etc.
                elif first_line.startswith(('HEADER', 'CRYST1', 'ATOM', 'HETATM', 'MODEL')):
                    return self.__class__.CONTENT_FLAG_PDB

            # Default to PDB if we can't determine
            return self.__class__.CONTENT_FLAG_PDB

        except Exception:
            return None

    def _introspect_content_flag(self):
        """
        Override base class to introspect PDB/mmCIF file using gemmi.

        Returns:
            int: CONTENT_FLAG_PDB (1), CONTENT_FLAG_MMCIF (2), or None if undetermined
        """
        import gemmi
        from pathlib import Path

        input_path = self.getFullPath()
        if not input_path or not Path(input_path).exists():
            return None

        try:
            # Use gemmi to read the structure - it auto-detects format
            structure = gemmi.read_structure(input_path)

            # Check the file extension to determine format
            path = Path(input_path)
            suffix = path.suffix.lower()

            # mmCIF files
            if suffix in ['.cif', '.mmcif']:
                return self.CONTENT_FLAG_MMCIF
            # PDB files
            elif suffix in ['.pdb', '.ent']:
                return self.CONTENT_FLAG_PDB
            else:
                # Default to PDB if ambiguous
                return self.CONTENT_FLAG_PDB

        except Exception as e:
            # If gemmi fails, fall back to simple text introspection
            try:
                with open(input_path, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith('data_') or first_line.startswith('loop_'):
                        return self.CONTENT_FLAG_MMCIF
                    else:
                        return self.CONTENT_FLAG_PDB
            except Exception:
                return None

    def fileExtensions(self):
        """
        Return appropriate file extension(s) for CPdbDataFile.

        CPdbDataFile can be either PDB or mmCIF format. We determine which by:
        1. If file exists: introspect using gemmi
        2. If contentFlag is set: use that (2 = mmCIF, other = PDB)
        3. Default: PDB

        Returns:
            list: [primary_extension] - either ['pdb'] or ['mmcif']
        """
        from pathlib import Path

        # Check if file exists and introspect it
        full_path = self.getFullPath()
        if full_path and Path(full_path).exists():
            try:
                # Use gemmi to determine format from actual file
                import gemmi
                structure = gemmi.read_structure(full_path)

                # Check file extension to determine format
                suffix = Path(full_path).suffix.lower()
                if suffix in ['.cif', '.mmcif']:
                    return ['mmcif']
                else:
                    return ['pdb']
            except Exception:
                # If gemmi fails, default to pdb
                return ['pdb']

        # For new files (not yet created), check contentFlag if set
        content_flag = 0
        if hasattr(self, 'contentFlag') and self.contentFlag is not None:
            if hasattr(self.contentFlag, 'value'):
                content_flag = self.contentFlag.value if self.contentFlag.value is not None else 0
            else:
                content_flag = int(self.contentFlag) if self.contentFlag else 0

        # contentFlag: 2 = mmCIF, anything else = PDB (default)
        # Note: mmCIF files use .cif extension (not .mmcif)
        if content_flag == self.CONTENT_FLAG_MMCIF:
            return ['cif']
        else:
            return ['pdb']

    def isMMCIF(self) -> bool:
        """
        Check if this coordinate file is in mmCIF format.

        This method loads the file content if necessary and checks the format.
        It provides a convenient API for plugins to detect mmCIF files.

        Returns:
            bool: True if file is mmCIF format, False otherwise

        Example:
            >>> pdb_file = CPdbDataFile()
            >>> pdb_file.setFullPath('/path/to/model.cif')
            >>> if pdb_file.isMMCIF():
            ...     print("File is in mmCIF format")
        """
        # First, try quick check based on file extension
        # Note: We can trust .cif/.mmcif extensions, but NOT .pdb extensions
        # (legacy wrappers may copy mmCIF files to .pdb paths)
        full_path = self.getFullPath()
        if full_path:
            from pathlib import Path
            suffix = Path(full_path).suffix.lower()
            if suffix in ['.cif', '.mmcif']:
                return True
            # Don't quick-fail on .pdb - it might actually be mmCIF content
            # Fall through to content checking below

        # Quick file content check (before trying to load with gemmi)
        # This is faster and avoids triggering loadFile() which might fail
        if full_path:
            from pathlib import Path
            if Path(full_path).exists():
                try:
                    # Quick check: read first line to see if it starts with mmCIF markers
                    with open(full_path, 'r') as f:
                        first_line = f.readline().strip()
                        # mmCIF files start with "data_" or have "loop_" structures early on
                        if first_line.startswith('data_'):
                            return True
                        # Also check second line for loop_ (some mmCIF files have comments first)
                        second_line = f.readline().strip()
                        if second_line.startswith('data_') or first_line.startswith('loop_') or second_line.startswith('loop_'):
                            return True
                except Exception:
                    pass  # If we can't read, fall through to other checks

        # If file content is already loaded (from previous loadFile call), check it
        # Note: We check this AFTER the quick file peek to avoid triggering loadFile()
        if hasattr(self, '_fileContent') and self._fileContent is not None:
            if hasattr(self._fileContent, 'isMMCIF'):
                return self._fileContent.isMMCIF()

        # Fallback: check contentFlag
        if hasattr(self, 'contentFlag') and self.contentFlag is not None:
            if hasattr(self.contentFlag, 'value'):
                content_flag = self.contentFlag.value if self.contentFlag.value is not None else 0
            else:
                content_flag = int(self.contentFlag) if self.contentFlag else 0
            return content_flag == self.CONTENT_FLAG_MMCIF

        return False

    def getSelectedAtomsPdbFile(self, fileName=None):
        """
        Apply atom selection and write selected atoms to file.

        This method replicates the legacy mmdb-based selection behavior,
        using our modern gemmi-based selection system.

        Args:
            fileName: Output file path (required)

        Returns:
            0 on success, non-zero on failure

        Example:
            >>> pdb_file = CPdbDataFile()
            >>> pdb_file.setFullPath('/path/to/model.pdb')
            >>> pdb_file.selection.text.set("A/27-50")
            >>> pdb_file.getSelectedAtomsPdbFile("/tmp/selected.pdb")
        """
        import shutil
        from pathlib import Path

        if fileName is None:
            raise ValueError("fileName parameter is required")

        # If no selection is set, just copy the file
        if not self.isSelectionSet():
            input_path = self.getFullPath()
            if input_path and Path(input_path).exists():
                shutil.copyfile(input_path, fileName)
                return 0
            else:
                print(f"Error: Input file does not exist: {input_path}")
                return 1

        # Load the file if not already loaded
        self.loadFile()

        # Get the selection string
        selection_string = str(self.selection.text.value) if hasattr(self.selection.text, 'value') else str(self.selection.text)

        try:
            # Interpret the selection
            n_atoms, selected_atoms = self.fileContent.interpretSelection(selection_string)

            # Write the selection to file
            rc = self.fileContent.writeSelection(selected_atoms, fileName)

            return rc

        except Exception as e:
            print(f"Error applying selection '{selection_string}': {e}")
            import traceback
            traceback.print_exc()
            return 1

    def convertFormat(self, toFormat, fileName):
        """
        Convert PDB/mmCIF file to specified format and write to fileName.

        Legacy compatibility method for wrappers that call convertFormat.
        Uses gemmi for format conversion.

        Args:
            toFormat: Target format ('pdb', 'cif', 'mmcif', etc.)
            fileName: Output file path (may be string, CString, or CFilePath)

        Raises:
            Exception: If file doesn't exist, can't be loaded, or can't be written
        """
        from ccp4i2.core.base_object.error_reporting import CException
        import os

        # Convert fileName to string if it's a CData object (CString, CFilePath, etc.)
        output_path = str(fileName)

        # Check if source file exists - try getFullPath() first, then str(self)
        source_path = self.getFullPath()
        if not source_path:
            source_path = str(self)

        if not source_path or not os.path.exists(source_path):
            raise CException(self.__class__, 410, f"Source file does not exist: {source_path}")

        # Load the file content if not already loaded
        self.loadFile()

        # Check if gemmi structure was loaded
        if not hasattr(self.fileContent, '_gemmi_structure') or self.fileContent._gemmi_structure is None:
            raise CException(self.__class__, 411, f"Failed to load structure from {source_path}")

        # Remove existing output file if present
        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except Exception as e:
                raise CException(self.__class__, 412, f"Cannot remove existing file {output_path}: {str(e)}")

        # Write in requested format using gemmi
        try:
            import gemmi

            if toFormat.lower() in ['cif', 'mmcif']:
                # Write as mmCIF
                self.fileContent._gemmi_structure.make_mmcif_document().write_file(output_path)
            else:
                # Write as PDB (default)
                pdb_str = self.fileContent._gemmi_structure.make_pdb_string()
                with open(output_path, 'w') as f:
                    f.write(pdb_str)

        except Exception as e:
            raise CException(self.__class__, 413, f"Failed to write {output_path}: {str(e)}")

        if not os.path.exists(output_path):
            raise CException(self.__class__, 413, f"Output file was not created: {output_path}")

    def replaceMSE(self, fileName):
        """
        Replace MSE (selenomethionine) residues with MET (methionine) and write to fileName.

        Legacy compatibility method for wrappers that call replaceMSE.
        Uses gemmi to load, modify, and write the structure.
        Preserves the input file format (PDB or mmCIF).

        This converts:
        - MSE residue name -> MET
        - SE atom name -> SD (sulfur delta)
        - HETATM records -> ATOM records

        Args:
            fileName: Output file path (may be string, CString, or CFilePath)

        Raises:
            Exception: If file doesn't exist, can't be loaded, or can't be written
        """
        from ccp4i2.core.base_object.error_reporting import CException
        import os

        # Convert fileName to string if it's a CData object (CString, CFilePath, etc.)
        output_path = str(fileName)

        # Check if source file exists - try getFullPath() first, then str(self)
        source_path = self.getFullPath()
        if not source_path:
            source_path = str(self)

        if not source_path or not os.path.exists(source_path):
            raise CException(self.__class__, 410, f"Source file does not exist: {source_path}")

        # Detect input format and load the structure using gemmi
        try:
            import gemmi

            # Try to detect if input is mmCIF or PDB
            is_mmcif = False
            try:
                gemmi.cif.read(source_path)
                is_mmcif = True
            except (ValueError, RuntimeError):
                # Not mmCIF, assume PDB
                is_mmcif = False

            # Load the structure
            structure = gemmi.read_structure(source_path)
        except Exception as e:
            raise CException(self.__class__, 411, f"Failed to load structure from {source_path}: {str(e)}")

        # Replace MSE with MET
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.name == 'MSE':
                        residue.name = 'MET'
                        # Rename SE atom to SD
                        for atom in residue:
                            if atom.name == 'SE' or atom.name.strip() == 'SE':
                                atom.name = 'SD'
                                atom.element = gemmi.Element('S')

        # Remove existing output file if present
        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except Exception as e:
                raise CException(self.__class__, 412, f"Cannot remove existing file {output_path}: {str(e)}")

        # Write the modified structure preserving the input format
        try:
            if is_mmcif:
                # Write as mmCIF
                structure.make_mmcif_document().write_file(output_path)
            else:
                # Write as PDB
                pdb_str = structure.make_pdb_string()
                with open(output_path, 'w') as f:
                    f.write(pdb_str)
        except Exception as e:
            raise CException(self.__class__, 413, f"Failed to write {output_path}: {str(e)}")

        if not os.path.exists(output_path):
            raise CException(self.__class__, 413, f"Output file was not created: {output_path}")


class CPdbDataFileList(CPdbDataFileListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CPdbDataFileListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPdbEnsembleItem(CPdbEnsembleItemStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    Extends CPdbEnsembleItemStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def isSet(self, attributeName=None, allowDefault=True):
        """
        Override isSet() to ensure ensemble items without a structure file
        are not considered "set" during XML serialization.

        An ensemble item is only meaningful if it has a structure file.
        This prevents partial/empty items from being saved to XML.
        """
        if attributeName is None:
            # Check if the entire object is set - requires structure to be set
            return hasattr(self, 'structure') and self.structure.isSet()
        else:
            # Check specific attribute using parent's implementation
            return super().isSet(attributeName, allowDefault=allowDefault)

    def getEtree(self, name: str = None, excludeUnset: bool = False, allSet: bool = False):
        """
        Override getEtree() to conditionally serialize identity_to_target OR rms_to_target.

        Legacy expectation: Only ONE of identity_to_target or rms_to_target should be
        written to XML, with preference given to whichever is non-default/non-zero.

        Args:
            name: Optional name for the XML element
            excludeUnset: If True, only serialize fields that have been explicitly set
            allSet: If True, only serialize if ALL registered attributes are set

        Returns:
            xml.etree.ElementTree.Element representing this object
        """
        import xml.etree.ElementTree as ET

        # Call parent to get the base element
        elem = super().getEtree(name=name, excludeUnset=excludeUnset, allSet=allSet)

        # Determine which similarity metric to keep (if any)
        identity_val = 0.0
        rms_val = 0.0
        identity_is_set = False
        rms_is_set = False

        if hasattr(self, 'identity_to_target') and self.identity_to_target is not None:
            identity_val = getattr(self.identity_to_target, 'value', 0.0)
            identity_is_set = self.identity_to_target.isSet()

        if hasattr(self, 'rms_to_target') and self.rms_to_target is not None:
            rms_val = getattr(self.rms_to_target, 'value', 0.0)
            rms_is_set = self.rms_to_target.isSet()

        # Decide which one to keep based on legacy rules:
        # - Keep identity if it's non-zero/non-default (preferred)
        # - Otherwise keep rms if it's non-zero/non-default
        # - Remove the other one from XML
        keep_identity = identity_is_set and identity_val != 0.0
        keep_rms = rms_is_set and rms_val != 0.0

        # If both are set and non-zero, prefer identity (legacy expectation)
        if keep_identity:
            # Remove rms_to_target from XML
            for child in elem.findall('rms_to_target'):
                elem.remove(child)
        elif keep_rms:
            # Remove identity_to_target from XML
            for child in elem.findall('identity_to_target'):
                elem.remove(child)
        # If neither is non-zero, keep both (or neither if excludeUnset=True handled by parent)

        return elem

    def validity(self):
        """
        Validate the ensemble item.

        Checks:
        - Structure file must be set
        - At least one of identity_to_target or rms_to_target must be set and non-zero

        Returns:
            CErrorReport containing validation errors/warnings
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport, SEVERITY_ERROR

        report = CErrorReport()

        # First, validate the structure file (call parent/child validity)
        if hasattr(self, 'structure') and self.structure is not None:
            structure_report = self.structure.validity()
            report.extend(structure_report)

        # Check that at least one similarity metric is set and non-zero
        identity_val = 0.0
        rms_val = 0.0

        if hasattr(self, 'identity_to_target') and self.identity_to_target is not None:
            identity_val = getattr(self.identity_to_target, 'value', 0.0) or 0.0

        if hasattr(self, 'rms_to_target') and self.rms_to_target is not None:
            rms_val = getattr(self.rms_to_target, 'value', 0.0) or 0.0

        has_valid_metric = (identity_val != 0.0) or (rms_val != 0.0)

        if not has_valid_metric:
            obj_path = self.object_path() if hasattr(self, 'object_path') else self.objectName()
            report.append(
                klass=self.__class__.__name__,
                code=100,
                details="Ensemble item requires either Identity or RMS to be set (non-zero)",
                name=obj_path,
                severity=SEVERITY_ERROR
            )

        return report


class CResidueRange(CResidueRangeStub):
    """
    A residue range selection
    
    Extends CResidueRangeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CResidueRangeList(CResidueRangeListStub):
    """
    A list of residue range selections
    
    Extends CResidueRangeListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSeqAlignDataFile(CSeqAlignDataFileStub):
    """
    A (multiple) sequence alignment file

    Extends CSeqAlignDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def identifyFile(self):
        """
        Identify the format of an alignment file using BioPython.

        Returns:
            tuple: (format_name, id_list) where format_name is the detected format
                   and id_list is list of sequence IDs in the alignment
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        import os

        err = CErrorReport()
        idList = []

        fileName = self.__str__()
        if not os.path.exists(fileName):
            err.append(self.__class__, 204, fileName, stack=False)
            return err, idList

        # Try common alignment formats with BioPython
        formats_to_try = ['clustal', 'fasta', 'phylip', 'stockholm']

        for fmt in formats_to_try:
            try:
                import Bio.AlignIO
                alignments = Bio.AlignIO.read(fileName, fmt)
                # Success - store format and extract IDs
                object.__setattr__(self, 'format', fmt)
                for record in alignments:
                    try:
                        idList.append(record.id)
                    except:
                        err.append(self.__class__, 205, fileName, stack=False)
                return err, idList
            except:
                continue

        # If no format worked, mark as unknown
        object.__setattr__(self, 'format', 'unknown')
        err.append(self.__class__, 204, fileName, stack=False)
        return err, idList

    def convertFormat(self, toFormat, fileName, reorder=None):
        """
        Convert alignment file to specified format and write to fileName.

        Legacy compatibility method for wrappers that call convertFormat.
        Uses BioPython for alignment format conversion.

        Args:
            toFormat: Target format ('clustal', 'fasta', 'phylip', 'stockholm', etc.)
            fileName: Output file path
            reorder: Optional reordering specification:
                    - None: Keep original order
                    - 'reverse': Reverse sequence order
                    - list: Reorder to specified indices

        Returns:
            CErrorReport with any errors encountered
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        import os

        # Identify source format if not already done
        if not hasattr(self, 'format') or self.format is None or self.format == 'unknown':
            self.identifyFile()

        if not hasattr(self, 'format') or self.format == 'unknown':
            return CErrorReport(self.__class__, 250, self.__str__() + ' to ' + str(fileName), stack=False)

        # Try reading the input file
        try:
            import Bio.AlignIO
            alignments = Bio.AlignIO.read(self.__str__(), self.format)
        except Exception as e:
            return CErrorReport(self.__class__, 202, f"{self.__str__()}: {str(e)}", stack=False)

        # Apply reordering if specified
        if reorder == 'reverse':
            from Bio.Align import MultipleSeqAlignment
            aliout = MultipleSeqAlignment([])
            for ii in range(len(alignments) - 1, -1, -1):
                aliout.append(alignments[ii])
            alignments = aliout
        elif isinstance(reorder, list):
            # Custom reordering by index list
            from Bio.Align import MultipleSeqAlignment
            aliout = MultipleSeqAlignment([])
            for idx in reorder:
                if 0 <= idx < len(alignments):
                    aliout.append(alignments[idx])
            alignments = aliout

        # Remove existing output file if present
        try:
            if os.path.exists(fileName):
                os.remove(fileName)
        except:
            return CErrorReport(self.__class__, 251, fileName)

        # Write output in requested format
        try:
            with open(fileName, "w") as out:
                Bio.AlignIO.write(alignments, out, toFormat)
        except Exception as e:
            return CErrorReport(self.__class__, 252, f"{fileName}: {str(e)}", stack=False)

        return CErrorReport()


class CSeqDataFile(CSeqDataFileStub):
    """
    A sequence file
    
    Extends CSeqDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSeqDataFileList(CSeqDataFileListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CSeqDataFileListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSequence(CSequenceStub, CBioPythonSeqInterface):
    """
    A string of sequence one-letter codes with BioPython loading support.

    This class represents a single biological sequence (protein or nucleic acid)
    with support for loading from various file formats using BioPython.

    Key features:
    - Load from FASTA, PIR, UniProt, and other formats
    - Automatic format detection
    - Molecular weight calculation using BioPython ProtParam
    - Save to FASTA format
    """

    def loadFile(self, fileName: str = None, format: str = 'unknown'):
        """
        Load sequence from file.

        Args:
            fileName: Path to sequence file. If None, gets path from parent CDataFile.
            format: Format hint ('internal', 'uniprot', 'fasta', 'pir', 'unknown')
        """
        from pathlib import Path
        from ccp4i2.core.base_object.error_reporting import CErrorReport

        # If no path provided, get from parent CDataFile
        if fileName is None:
            parent = self.get_parent()
            if parent is not None and hasattr(parent, 'getFullPath'):
                fileName = parent.getFullPath()

        # Return early if still no file path
        if not fileName:
            return

        if format == 'internal':
            # Internal format: simple FASTA-like with pipe-separated metadata
            self._loadInternalFile(fileName)
        elif format == 'uniprot':
            # UniProt XML format - not implemented in modern version
            err = CErrorReport()
            err.append(klass='CSequence', code=402,
                      details='UniProt XML format not yet supported in modern implementation')
            return err
        else:
            # Use BioPython interface
            data = self.loadExternalFile(fileName, format=format, record=0)
            if data is not None:
                # Populate CData attributes from loaded data
                if 'identifier' in data:
                    self.identifier = data['identifier']
                if 'sequence' in data:
                    self.sequence = data['sequence']
                    # Detect and set moleculeType based on sequence content
                    self.moleculeType = self._detect_molecule_type(data['sequence'])
                if 'name' in data:
                    self.name = data['name']
                if 'description' in data:
                    self.description = data['description']
                if 'referenceDb' in data:
                    self.referenceDb = data['referenceDb']
                if 'reference' in data:
                    self.reference = data['reference']

    def _detect_molecule_type(self, sequence: str) -> str:
        """
        Detect whether a sequence is PROTEIN, DNA, or RNA.

        Args:
            sequence: The sequence string (single-letter codes)

        Returns:
            "PROTEIN", "DNA", or "RNA"

        Notes:
            - Checks for presence of U (uracil) to identify RNA
            - Checks for T (thymine) to identify DNA
            - Defaults to PROTEIN for amino acid sequences
        """
        seq_upper = str(sequence).upper().replace('*', '').strip()

        # Remove common modifications and gaps
        seq_clean = seq_upper.replace('-', '').replace('.', '').replace(' ', '')

        if not seq_clean:
            return "PROTEIN"

        # Count nucleotide-specific characters
        has_u = 'U' in seq_clean
        has_t = 'T' in seq_clean

        # DNA/RNA detection
        dna_bases = set('ACGT')
        rna_bases = set('ACGU')

        # If majority are DNA/RNA bases, classify as nucleic acid
        base_count = sum(1 for c in seq_clean if c in dna_bases.union(rna_bases))
        base_fraction = base_count / len(seq_clean) if seq_clean else 0

        if base_fraction > 0.7:  # >70% nucleotide bases
            if has_u and not has_t:
                return "RNA"
            elif has_t and not has_u:
                return "DNA"
            elif has_u:  # Both U and T present - prefer RNA
                return "RNA"
            else:
                return "DNA"

        # Default to protein
        return "PROTEIN"

    def _loadInternalFile(self, fileName: str):
        """
        Load internal format sequence file.

        Internal format is simple FASTA with pipe-separated metadata:
        >referenceDb|reference|identifier
        sequence...
        """
        from pathlib import Path

        try:
            text = Path(fileName).read_text()
        except Exception as e:
            print(f'ERROR loading sequence file: {e}')
            return

        lines = text.split('\n')
        try:
            # Parse header line
            if len(lines[0]) > 0 and lines[0][0] == '>':
                header = lines[0][1:]
                splitList = header.split('|')

                if len(splitList) == 3:
                    self.referenceDb = splitList[0]
                    if len(splitList[1].strip()) > 0:
                        self.reference = splitList[1]
                    self.identifier = splitList[2]
                else:
                    self.identifier = header

                # Parse sequence (concatenate all remaining lines)
                seq = ''.join(lines[1:])
                self.sequence = seq
                # Detect and set moleculeType based on sequence content
                self.moleculeType = self._detect_molecule_type(seq)
        except Exception as e:
            print(f'ERROR parsing sequence file: {e}')

    def saveFile(self, fileName: str):
        """
        Save sequence to FASTA format file.

        Args:
            fileName: Output file path
        """
        from pathlib import Path

        # Build FASTA format
        text = '>' + str(self.identifier) + '\n' + str(self.sequence)

        # Write to file
        Path(fileName).write_text(text)

    def getAnalysis(self, mode: str = 'molecularWeight'):
        """
        Perform sequence analysis using BioPython.

        Args:
            mode: Analysis type ('molecularWeight' supported)

        Returns:
            Molecular weight in Daltons, or 0 if sequence not set or analysis fails
        """
        import re

        if mode == 'molecularWeight':
            if not self.sequence.isSet():
                return 0

            if not BIOPYTHON_AVAILABLE:
                print('BioPython not available for molecular weight calculation')
                return 0

            try:
                from Bio.SeqUtils.ProtParam import ProteinAnalysis

                # Remove non-standard amino acids for BioPython compatibility
                seq_str = str(self.sequence)
                seq_clean = re.sub('[^GALMFWKQESPVICYHRNDT]', '', seq_str)

                if len(seq_clean) == 0:
                    return 0

                pa = ProteinAnalysis(seq_clean)
                return pa.molecular_weight()
            except Exception as e:
                print(f'Error calculating molecular weight: {e}')
                return 0

        return 0

    def guiLabel(self) -> str:
        """
        Get display label for GUI.

        Returns:
            Identifier, or first 20 chars of sequence, or object name
        """
        if self.identifier.isSet():
            return str(self.identifier)
        elif self.sequence.isSet():
            seq_str = str(self.sequence)
            return seq_str[0:20] if len(seq_str) > 20 else seq_str
        else:
            return self.object_name() if hasattr(self, 'object_name') else 'CSequence'


class CSequenceAlignment(CSequenceAlignmentStub, CBioPythonSeqInterface):
    """
    An alignment of two or more sequences with BioPython AlignIO support.

    This class represents a multiple sequence alignment with support for
    loading from various alignment file formats using BioPython.

    Key features:
    - Load from CLUSTAL, FASTA, PHYLIP, Stockholm, Nexus, MSF formats
    - Automatic format detection
    - Handles both consecutive and interleaved alignment formats
    - Preserves gap characters for alignment information

    The alignment contains gaps ("-" characters) that are relevant to the alignment.
    Each aligned sequence is similar to CSequence but includes gap positions.
    """

    def loadFile(self, fileName: str, format: str = 'unknown'):
        """
        Load sequence alignment from file using BioPython AlignIO.

        Args:
            fileName: Path to alignment file
            format: Format hint ('clustal', 'fasta', 'phylip', 'stockholm', 'nexus', 'msf', 'unknown')
        """
        from pathlib import Path
        from ccp4i2.core.base_object.error_reporting import CErrorReport

        if not Path(fileName).exists():
            err = CErrorReport()
            err.append(klass='CSequenceAlignment', code=401,
                      details=f'Alignment file does not exist: {fileName}')
            return err

        if not BIOPYTHON_AVAILABLE:
            err = CErrorReport()
            err.append(klass='CSequenceAlignment', code=402,
                      details='BioPython not available - cannot load alignment file')
            return err

        # Use BioPython interface to load alignment
        data = self.loadExternalFile(fileName, format=format, record=0)

        if data is not None:
            # Populate CData attributes from loaded data
            if 'identifier' in data:
                self.identifier = data['identifier']

            # Store format information if available
            if 'format' in data:
                # Store in a private attribute for later use
                object.__setattr__(self, '_loaded_format', data['format'])

        return CErrorReport()  # Return empty error report on success

    def getSequenceCount(self) -> int:
        """
        Get the number of sequences in the alignment.

        Returns:
            Number of sequences, or 0 if not loaded
        """
        # This would require storing the full alignment data
        # For now, return 0 as a placeholder
        return 0

    def getAlignmentLength(self) -> int:
        """
        Get the length of the alignment (including gaps).

        Returns:
            Alignment length, or 0 if not loaded
        """
        # This would require storing the full alignment data
        # For now, return 0 as a placeholder
        return 0


class CSequenceMeta(CSequenceMetaStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CSequenceMetaStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSequenceString(CSequenceStringStub):
    """
    A string
    
    Extends CSequenceStringStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CTLSDataFile(CTLSDataFileStub):
    """
    A refmac TLS file
    
    Extends CTLSDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

