"""
Convert sequence files (FASTA, PIR, SEQ) to CAsuDataFile XML format.

This module provides utilities to read sequence files in various formats
and generate CCP4i2 ASU (Asymmetric Unit) content XML files that can be
used as input for crystallographic pipelines.

The ASU XML format contains sequence information with metadata:
- sequence: amino acid or nucleotide sequence
- nCopies: stoichiometry (number of copies in ASU)
- polymerType: PROTEIN, DNA, or RNA
- name: sequence identifier
- description: optional description

Example ASU XML structure:
    <ccp4i2_body>
      <seqList>
        <CAsuContentSeq>
          <sequence>MHHHHHHLVPR...</sequence>
          <nCopies>1</nCopies>
          <polymerType>PROTEIN</polymerType>
          <name>MyProtein</name>
          <description>My protein description</description>
        </CAsuContentSeq>
      </seqList>
    </ccp4i2_body>
"""

import logging
from pathlib import Path
from typing import Union, List, Dict, Optional
from datetime import datetime
import xml.etree.ElementTree as ET
import socket
import getpass

logger = logging.getLogger(__name__)


def detect_polymer_type(sequence: str) -> str:
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
    seq_upper = sequence.upper().replace('*', '').strip()

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


def parse_sequence_file(file_path: Union[str, Path]) -> List[Dict[str, str]]:
    """
    Parse a sequence file (FASTA, PIR, or other formats) using BioPython.

    Args:
        file_path: Path to sequence file

    Returns:
        List of dicts with keys: 'name', 'description', 'sequence', 'polymerType'

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file cannot be parsed or contains no sequences

    Example:
        >>> sequences = parse_sequence_file('gamma.pir')
        >>> sequences[0]['name']
        'gamma'
        >>> sequences[0]['polymerType']
        'PROTEIN'
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Sequence file not found: {file_path}")

    try:
        from Bio import SeqIO
    except ImportError:
        raise ImportError(
            "BioPython is required for sequence file parsing. "
            "Install with: pip install biopython"
        )

    sequences = []

    # Try different formats
    formats_to_try = ['fasta', 'pir', 'genbank', 'embl']

    parsed_records = None
    successful_format = None

    for fmt in formats_to_try:
        try:
            records = list(SeqIO.parse(str(file_path), fmt))
            if records:
                parsed_records = records
                successful_format = fmt
                logger.debug(f"Successfully parsed {file_path.name} as {fmt} format")
                break
        except Exception as e:
            logger.debug(f"Failed to parse as {fmt}: {e}")
            continue

    if not parsed_records:
        raise ValueError(
            f"Could not parse {file_path.name} as any supported format "
            f"({', '.join(formats_to_try)})"
        )

    logger.info(f"Parsed {len(parsed_records)} sequence(s) from {file_path.name} ({successful_format} format)")

    for record in parsed_records:
        # Get sequence string (remove terminal * if present)
        seq_str = str(record.seq).replace('*', '').strip()

        # Detect polymer type
        polymer_type = detect_polymer_type(seq_str)

        # Extract name and description
        name = record.id if record.id else "Unknown"
        description = record.description if record.description else name

        # If description starts with the ID, remove it to avoid redundancy
        if description.startswith(name):
            description = description[len(name):].strip()

        if not description:
            description = name

        sequences.append({
            'name': name,
            'description': description,
            'sequence': seq_str,
            'polymerType': polymer_type
        })

    return sequences


def create_asu_xml(
    sequences: List[Dict[str, str]],
    output_path: Optional[Union[str, Path]] = None,
    project_name: str = "i2run_project",
    project_id: str = "00000000000000000000000000000000"
) -> ET.Element:
    """
    Create a CAsuDataFile XML element from sequence information.

    Args:
        sequences: List of dicts with 'name', 'description', 'sequence', 'polymerType'
        output_path: Optional path to write XML file
        project_name: Project name for XML header
        project_id: Project ID for XML header (32-char hex)

    Returns:
        XML Element (root of the ASU XML tree)

    Example:
        >>> seqs = [{'name': 'gamma', 'description': 'gamma protein',
        ...          'sequence': 'MHHHHH...', 'polymerType': 'PROTEIN'}]
        >>> root = create_asu_xml(seqs, 'output.asu.xml')
    """
    # Create root element with namespace
    root = ET.Element('ccp4i2')
    root.set('xmlns:ccp4', 'http://www.ccp4.ac.uk/ccp4ns')

    # Create header
    header = ET.SubElement(root, 'ccp4i2_header')

    ET.SubElement(header, 'function').text = 'ASUCONTENT'

    try:
        user_id = getpass.getuser()
    except Exception:
        user_id = 'unknown'
    ET.SubElement(header, 'userId').text = user_id

    try:
        host_name = socket.gethostname()
    except Exception:
        host_name = 'localhost'
    ET.SubElement(header, 'hostName').text = host_name

    # Format: "HH:MM DD/Mon/YY"
    creation_time = datetime.now().strftime('%H:%M %d/%b/%y')
    ET.SubElement(header, 'creationTime').text = creation_time

    ET.SubElement(header, 'projectName').text = project_name
    ET.SubElement(header, 'projectId').text = project_id

    # Create body
    body = ET.SubElement(root, 'ccp4i2_body')
    seq_list = ET.SubElement(body, 'seqList')

    # Add each sequence
    for seq_info in sequences:
        seq_elem = ET.SubElement(seq_list, 'CAsuContentSeq')

        ET.SubElement(seq_elem, 'sequence').text = seq_info['sequence']
        ET.SubElement(seq_elem, 'nCopies').text = '1'
        ET.SubElement(seq_elem, 'polymerType').text = seq_info['polymerType']
        ET.SubElement(seq_elem, 'name').text = seq_info['name']
        ET.SubElement(seq_elem, 'description').text = seq_info['description']

    # Write to file if requested
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Pretty print with indentation
        _indent(root)

        tree = ET.ElementTree(root)
        tree.write(
            str(output_path),
            encoding='ASCII',
            xml_declaration=True
        )
        logger.info(f"Created ASU XML file: {output_path}")

    return root


def _indent(elem, level=0):
    """
    Add indentation to XML element for pretty printing.

    This is a compatibility function for Python < 3.9 which doesn't have ET.indent().
    """
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for child in elem:
            _indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def convert_sequence_file_to_asu(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    project_name: str = "i2run_project",
    project_id: str = "00000000000000000000000000000000"
) -> Path:
    """
    Convert a sequence file to CAsuDataFile XML format.

    This is the main convenience function that combines parsing and XML generation.

    Args:
        input_file: Path to input sequence file (FASTA, PIR, etc.)
        output_file: Path for output XML file (default: input_file.asu.xml)
        project_name: Project name for XML header
        project_id: Project ID for XML header

    Returns:
        Path to generated XML file

    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If file cannot be parsed

    Example:
        >>> asu_path = convert_sequence_file_to_asu('gamma.pir')
        >>> print(asu_path)
        gamma.asu.xml
    """
    input_file = Path(input_file)

    if output_file is None:
        # Default: same directory, same name with .asu.xml extension
        output_file = input_file.with_suffix('.asu.xml')
    else:
        output_file = Path(output_file)

    # Parse sequences from input file
    sequences = parse_sequence_file(input_file)

    if not sequences:
        raise ValueError(f"No sequences found in {input_file}")

    # Create ASU XML
    create_asu_xml(
        sequences=sequences,
        output_path=output_file,
        project_name=project_name,
        project_id=project_id
    )

    return output_file
