"""
CCP4 Utility Functions

Collection of utility functions for CCP4 crystallographic operations.
This module has no CData dependencies - pure utility functions.
"""

from pathlib import Path
from typing import List, Dict, Union, Optional
import gemmi
import numpy as np


class MtzMergeError(Exception):
    """Errors during MTZ merging operations."""
    pass


class MtzSplitError(Exception):
    """Errors during MTZ splitting operations."""
    pass


def split_mtz_file(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    column_mapping: dict
) -> Path:
    """
    Split/extract columns from an MTZ file using gemmi (completely CData-agnostic).

    Creates a new MTZ file with a complete unique reflection set for the space group,
    containing only the specified columns. Uses gemmi.make_miller_array() to ensure
    all expected reflections (including systematic absences) are present.

    Args:
        input_path: Path to input MTZ file
        output_path: Path for output MTZ file
        column_mapping: Dictionary mapping input column names to output column names
                       e.g., {'FMEAN': 'F', 'SIGFMEAN': 'SIGF'}

    Returns:
        Path: Full path to created output file

    Raises:
        FileNotFoundError: If input MTZ file doesn't exist
        ValueError: If column not found in input
        MtzSplitError: If gemmi operations fail

    Example:
        >>> split_mtz_file(
        ...     input_path='/data/full.mtz',
        ...     output_path='/data/mini.mtz',
        ...     column_mapping={'FMEAN': 'F', 'SIGFMEAN': 'SIGF'}
        ... )
        Path('/data/mini.mtz')
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    # Validate input
    if not input_path.exists():
        raise FileNotFoundError(f"Input MTZ file not found: {input_path}")

    if not column_mapping:
        raise ValueError("column_mapping cannot be empty")

    try:
        # Read input MTZ and ensure reflections are in ASU
        mtzin = gemmi.read_mtz_file(str(input_path))
        mtzin.ensure_asu()

        # Validate all requested columns exist
        available_columns = mtzin.column_labels()
        for input_col in column_mapping.keys():
            if input_col not in available_columns:
                raise ValueError(
                    f"Column '{input_col}' not found in {input_path}. "
                    f"Available columns: {available_columns}"
                )

        # Create output MTZ with complete unique reflection set
        mtzout = gemmi.Mtz()
        mtzout.spacegroup = mtzin.spacegroup
        mtzout.cell = mtzin.cell

        # Add HKL_base dataset for H, K, L columns
        hkl_base = mtzout.add_dataset('HKL_base')
        mtzout.add_column('H', 'H')
        mtzout.add_column('K', 'H')
        mtzout.add_column('L', 'H')

        # Create complete unique reflection set for the space group
        # This ensures all expected reflections (including absences) are present
        uniques = gemmi.make_miller_array(
            mtzout.cell,
            mtzout.spacegroup,
            mtzin.resolution_high(),
            mtzin.resolution_low()
        )
        mtzout.set_data(uniques)

        # Determine if we need a data dataset (for non-HKL columns)
        dataset = hkl_base
        for input_col_name in column_mapping.keys():
            col = mtzin.column_with_label(input_col_name)
            if col.dataset_id > 0 and len(mtzin.datasets) > 1:
                src_dataset = mtzin.dataset(col.dataset_id)
                dataset = mtzout.add_dataset(src_dataset.project_name)
                dataset.crystal_name = src_dataset.crystal_name
                dataset.dataset_name = src_dataset.dataset_name
                dataset.wavelength = src_dataset.wavelength
                break

        # Copy data columns using copy_column for automatic H,K,L matching
        for input_col_name, output_col_name in column_mapping.items():
            src_col = mtzin.column_with_label(input_col_name)
            # Use copy_column(-1, col) to automatically match H,K,L indices
            new_col = mtzout.copy_column(-1, src_col)
            new_col.label = output_col_name
            new_col.dataset_id = dataset.id

        # Set history
        mtzout.history = [
            f'MTZ file created from {input_path.name} using split_mtz_file (gemmi)',
            f'Columns: {", ".join(f"{i}->{o}" for i, o in column_mapping.items())}'
        ]

        # Write output file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        mtzout.write_to_file(str(output_path))

        if not output_path.exists():
            raise MtzSplitError(f"Output file not created: {output_path}")

        return output_path

    except (ValueError, FileNotFoundError):
        raise
    except Exception as e:
        raise MtzSplitError(f"Failed to split MTZ file: {e}")


def merge_mtz_files_cad(
    input_specs: List[dict],
    output_path: Union[str, Path],
    merge_strategy: str = 'first'
) -> Path:
    """
    Merge multiple MTZ files using CAD (CCP4's official merging tool).

    This function uses the external `cad` command from CCP4 to properly merge
    MTZ files with different reflection sets, handling H,K,L matching correctly.

    Args: Same as merge_mtz_files
    Returns: Path to merged MTZ file
    """
    import subprocess
    import os
    import shutil

    output_path = Path(output_path)

    # Check if CAD is available
    cad_exe = shutil.which('cad')
    if not cad_exe:
        cbin = os.environ.get('CBIN')
        if cbin:
            cad_path = os.path.join(cbin, 'cad')
            if os.path.exists(cad_path):
                cad_exe = cad_path

    if not cad_exe:
        raise MtzMergeError("CAD executable not found. Please ensure CCP4 is installed and sourced.")

    # Build CAD command
    cmd = [cad_exe, 'HKLOUT', str(output_path)]

    # Add input files
    for i, spec in enumerate(input_specs, 1):
        cmd.extend(['HKLIN' + str(i), str(spec['path'])])

    # Build labin commands for column selection/renaming
    labin_lines = []
    for i, spec in enumerate(input_specs, 1):
        column_mapping = spec.get('column_mapping', {})
        if column_mapping:
            parts = [f'E{i}={in_col}' for in_col, out_col in column_mapping.items()]
            labin_lines.append(f'labin file {i} {" ".join(parts)}')

    # Run CAD
    stdin_text = '\n'.join(labin_lines + ['END'])

    try:
        result = subprocess.run(
            cmd,
            input=stdin_text,
            capture_output=True,
            text=True,
            timeout=60,
            env=os.environ.copy()
        )

        if result.returncode != 0:
            raise MtzMergeError(f"CAD failed: {result.stderr}")

        if not output_path.exists():
            raise MtzMergeError(f"CAD did not create output file: {output_path}")

        return output_path

    except subprocess.TimeoutExpired:
        raise MtzMergeError("CAD timed out")
    except Exception as e:
        raise MtzMergeError(f"CAD execution failed: {e}")


def merge_mtz_files(
    input_specs: List[dict],
    output_path: Union[str, Path],
    merge_strategy: str = 'first'
) -> Path:
    """
    Merge multiple MTZ files using CAD from CCP4 (CData-agnostic).

    This is a low-level utility that merges reflection data from multiple
    MTZ files into a single output file. It has NO knowledge of CMiniMtzDataFile,
    CONTENT_SIGNATURE_LIST, or any CData conventions.

    Args:
        input_specs: List of dictionaries, each containing:
            {
                'path': str/Path,                    # Filesystem path to input MTZ
                'column_mapping': Dict[str, str]     # input_label -> output_label
            }

            The column_mapping dictionary specifies which columns to copy and what
            to call them in the output. For example:
                {'F': 'F_NAT', 'SIGF': 'SIGF_NAT'}  # Copy F as F_NAT, SIGF as SIGF_NAT
                {'FreeR_flag': 'FreeR_flag'}          # Copy FreeR_flag unchanged

        output_path: Where to write the merged MTZ file

        merge_strategy: How to handle output column name conflicts:
            - 'first': Keep column from first file (default)
            - 'last': Keep column from last file (not fully implemented)
            - 'error': Raise error on conflicts
            - 'rename': Auto-rename conflicts (F, F_1, F_2, ...)

    Returns:
        Path: Full path to created output file

    Raises:
        FileNotFoundError: If input MTZ file doesn't exist
        ValueError: If column not found or conflict with strategy='error'
        MtzMergeError: If gemmi operations fail

    Example:
        >>> merge_mtz_files(
        ...     input_specs=[
        ...         {
        ...             'path': '/data/native.mtz',
        ...             'column_mapping': {'F': 'F_NAT', 'SIGF': 'SIGF_NAT'}
        ...         },
        ...         {
        ...             'path': '/data/free.mtz',
        ...             'column_mapping': {'FreeR_flag': 'FreeR_flag'}
        ...         }
        ...     ],
        ...     output_path='/data/merged.mtz',
        ...     merge_strategy='first'
        ... )
        Path('/data/merged.mtz')
    """
    output_path = Path(output_path)

    # Validate inputs
    if not input_specs:
        raise ValueError("input_specs cannot be empty")

    # Read first file to get base metadata and reflections
    first_spec = input_specs[0]
    first_path = Path(first_spec['path'])

    if not first_path.exists():
        raise FileNotFoundError(f"Input MTZ file not found: {first_path}")

    try:
        first_mtz = gemmi.read_mtz_file(str(first_path))
    except Exception as e:
        raise MtzMergeError(f"Failed to read {first_path}: {e}")

    # Ensure first MTZ has reflections in ASU
    first_mtz.ensure_asu()

    # Create new output MTZ with complete reflection set for the space group
    out_mtz = gemmi.Mtz()
    out_mtz.spacegroup = first_mtz.spacegroup
    out_mtz.cell = first_mtz.cell

    # Add HKL_base dataset for H, K, L columns
    hkl_base = out_mtz.add_dataset('HKL_base')
    out_mtz.add_column('H', 'H')
    out_mtz.add_column('K', 'H')
    out_mtz.add_column('L', 'H')

    # Create complete unique reflection set for the space group
    # This ensures all files will have a common reflection list
    # Note: resolution_high() returns HIGH resolution (small d-spacing)
    #       resolution_low() returns LOW resolution (large d-spacing)
    #       make_miller_array expects: (cell, spacegroup, d_min, d_max)
    #       So d_min should be resolution_high() and d_max should be resolution_low()
    uniques = gemmi.make_miller_array(
        out_mtz.cell,
        out_mtz.spacegroup,
        first_mtz.resolution_high(),  # d_min (high resolution, small value)
        first_mtz.resolution_low()     # d_max (low resolution, large value)
    )
    out_mtz.set_data(uniques)

    # Track which columns have been added to detect conflicts
    added_columns = set()

    # Process each input file
    for spec_idx, spec in enumerate(input_specs):
        input_path = Path(spec['path'])
        column_mapping = spec.get('column_mapping', {})

        # Validate required keys
        if 'path' not in spec or 'column_mapping' not in spec:
            raise ValueError(f"input_specs[{spec_idx}] missing required key 'path' or 'column_mapping'")

        if not column_mapping:
            continue  # Skip empty column mapping

        # Check file exists
        if not input_path.exists():
            raise FileNotFoundError(f"Input MTZ file not found: {input_path}")

        # Read MTZ file
        try:
            in_mtz = gemmi.read_mtz_file(str(input_path))
            in_mtz.ensure_asu()  # Ensure reflections are in ASU
        except Exception as e:
            raise MtzMergeError(f"Failed to read {input_path}: {e}")

        # Validate compatibility (same space group and cell)
        if in_mtz.spacegroup.number != out_mtz.spacegroup.number:
            raise MtzMergeError(
                f"Incompatible space groups: {first_path} has {out_mtz.spacegroup.hm}, "
                f"{input_path} has {in_mtz.spacegroup.hm}"
            )

        # Check cell parameters (allow small differences)
        cell_diff = sum(abs(a - b) for a, b in zip(in_mtz.cell.parameters, out_mtz.cell.parameters))
        if cell_diff > 0.1:  # Tolerance for cell parameter differences
            raise MtzMergeError(
                f"Incompatible unit cells: {first_path} has {out_mtz.cell.parameters}, "
                f"{input_path} has {in_mtz.cell.parameters}"
            )

        # Add dataset if needed (for data columns, not H,K,L)
        if len(out_mtz.datasets) < 2:
            # Get first data dataset from source
            src_dataset = in_mtz.datasets[1] if len(in_mtz.datasets) > 1 else in_mtz.datasets[0]
            ds = out_mtz.add_dataset(src_dataset.dataset_name or 'data')
            ds.crystal_name = src_dataset.crystal_name
            ds.wavelength = src_dataset.wavelength

        # Copy requested columns (input_label -> output_label)
        for input_label, output_label in column_mapping.items():
            # Find column in input MTZ
            src_col = in_mtz.column_with_label(input_label)
            if src_col is None:
                raise ValueError(
                    f"Column '{input_label}' not found in {input_path}. "
                    f"Available columns: {[c.label for c in in_mtz.columns]}"
                )

            # Check for conflicts
            if output_label in added_columns:
                if merge_strategy == 'error':
                    raise ValueError(
                        f"Column conflict: '{output_label}' already exists in output. "
                        f"Source: {input_path}, column: {input_label}"
                    )
                elif merge_strategy == 'first':
                    # Skip this column, keep the first one
                    continue
                elif merge_strategy == 'last':
                    # Remove existing column and add new one
                    # Gemmi doesn't support removing columns, so we'll skip for now
                    #  TODO: implement column replacement
                    pass
                elif merge_strategy == 'rename':
                    # Auto-rename with suffix
                    counter = 1
                    new_label = f"{output_label}_{counter}"
                    while new_label in added_columns:
                        counter += 1
                        new_label = f"{output_label}_{counter}"
                    output_label = new_label

            # Copy column using gemmi's automatic H,K,L matching
            try:
                # Use copy_column with -1 to automatically match H,K,L indices
                new_col = out_mtz.copy_column(-1, src_col)
                new_col.label = output_label
                new_col.dataset_id = 1  # Use dataset 1 for all data columns

                added_columns.add(output_label)

            except Exception as e:
                raise MtzMergeError(f"Failed to copy column '{input_label}' from {input_path}: {e}")

    # Update resolution limits
    out_mtz.update_reso()

    # Add history
    history_lines = [
        'Merged MTZ files using merge_mtz_files (gemmi)',
        f'Input files: {len(input_specs)}',
        f'Strategy: {merge_strategy}'
    ]
    for spec in input_specs:
        mapping = spec.get('column_mapping', {})
        cols = ', '.join(f"{in_col}->{out_col}" for in_col, out_col in mapping.items())
        history_lines.append(f"  - {spec['path']}: {cols}")
    out_mtz.history = history_lines

    # Write output file
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_mtz.write_to_file(str(output_path))

        # Explicitly sync file to disk to ensure MTZ history is fully written
        # This prevents UnicodeDecodeError when servalcat/gemmi reads MTZ immediately after
        import os
        with open(output_path, 'rb') as f:
            os.fsync(f.fileno())
    except Exception as e:
        raise MtzMergeError(f"Failed to write output MTZ to {output_path}: {e}")

    return output_path


def getCCP4Dir():
    """
    Get the CCP4 installation directory.

    Returns the root directory of the CCP4 suite installation from the $CCP4
    environment variable. This is used by plugins to locate CCP4 resources
    like monomer libraries, reference data, etc.

    Returns:
        str: Absolute path to CCP4 root directory, or '.' if not set

    Example:
        >>> ccp4_root = getCCP4Dir()
        >>> # Build path to monomer library
        >>> comp_id = 'ATP'
        >>> cif_path = os.path.join(ccp4_root, 'lib', 'data', 'monomers',
        ...                         comp_id[0].lower(), comp_id + '.cif')

    Notes:
        - Reads from $CCP4 environment variable
        - Returns '.' (current directory) if $CCP4 is not set (due to normpath behavior)
        - Normalizes the path before returning
        - Matches legacy ccp4i2 behavior exactly
    """
    import os

    try:
        target = os.path.normpath(os.environ.get('CCP4', ''))
    except Exception:
        target = ''
    return target


def getTMP(**kw):
    """
    Get the temporary directory for CCP4 jobs.

    Returns:
        str: Path to temporary directory

    Notes:
        - Legacy compatibility function for plugins that call getTMP()
        - Returns value from $TMP, $TEMP, or $TMPDIR environment variables
        - Falls back to '/tmp' on Unix or 'C:\\Temp' on Windows
        - Creates the directory if it doesn't exist
    """
    import os
    import tempfile

    # Try standard environment variables in order of preference
    tmp_dir = os.environ.get('TMP') or os.environ.get('TEMP') or os.environ.get('TMPDIR')

    if not tmp_dir:
        # Fall back to system default
        tmp_dir = tempfile.gettempdir()

    # Ensure directory exists
    if tmp_dir and not os.path.exists(tmp_dir):
        try:
            os.makedirs(tmp_dir, exist_ok=True)
        except Exception:
            pass

    return os.path.normpath(tmp_dir) if tmp_dir else tempfile.gettempdir()


def getCCP4I2Dir(**kw):
    """
    Get the CCP4i2 installation directory.

    Returns the root directory of the CCP4i2 installation by examining
    the location of this module and going up one directory level.

    This is used by plugins to locate resources like smartie, reports, etc.

    Returns:
        str: Absolute path to CCP4i2 root directory

    Example:
        >>> ccp4i2_root = getCCP4I2Dir()
        >>> smartie_path = os.path.join(ccp4i2_root, 'smartie')
    """
    import os
    import sys

    # Use the module's file location to find CCP4i2 root
    # CCP4Utils is at <CCP4I2_ROOT>/core/CCP4Utils.py
    # So we go up one level from core/ to get CCP4I2_ROOT
    f = os.path.normpath(__import__('core.CCP4Utils').__file__)
    return os.path.split(os.path.split(f)[0])[0]


def writeXML(file, xml_string):
    """
    Write XML string to file.

    Legacy compatibility function used by plugins to write XML output.

    Args:
        file: File object opened for writing
        xml_string: XML content as string or bytes
    """
    if isinstance(xml_string, bytes):
        xml_string = xml_string.decode('utf-8')
    file.write(xml_string)


def openFileToEtree(fileName: Union[str, Path], printout: bool = False, useLXML: bool = True):
    """
    Parse an XML file and return an ElementTree element.

    Legacy compatibility function used by plugins to read XML files.
    Reads the file content and parses it to avoid segfaults that can occur
    with direct etree.parse() on some systems.

    Args:
        fileName: Path to XML file to read
        printout: If True, print the parsed XML tree
        useLXML: If True, use lxml parser; if False, use xml.etree.ElementTree

    Returns:
        ElementTree element representing the parsed XML

    Raises:
        Exception: If file cannot be read or parsed
    """
    import os

    # Normalize path and read file
    file_path = os.path.normpath(str(fileName))
    with open(file_path, 'r', encoding='utf-8') as f:
        xml_content = f.read()

    # Parse using appropriate parser
    if useLXML:
        try:
            from lxml import etree
            xml_bytes = xml_content.encode('utf-8')
            tree = etree.fromstring(xml_bytes)
        except ImportError:
            # Fall back to standard library if lxml not available
            import xml.etree.ElementTree as ET
            tree = ET.fromstring(xml_content)
    else:
        import xml.etree.ElementTree as ET
        tree = ET.fromstring(xml_content)

    if printout:
        try:
            from lxml import etree
            print(etree.tostring(tree, pretty_print=True, encoding='unicode'))
        except ImportError:
            import xml.etree.ElementTree as ET
            print(ET.tostring(tree, encoding='unicode'))

    return tree


def readFile(fileName: Union[str, Path]) -> str:
    """
    Read a text file and return its contents as a string.

    Legacy compatibility function used by plugins to read text files.

    Args:
        fileName: Path to the file to read

    Returns:
        File contents as string

    Raises:
        FileNotFoundError: If file doesn't exist
        IOError: If file cannot be read

    Example:
        >>> content = readFile('/path/to/file.txt')
        >>> lines = content.split('\n')
    """
    file_path = Path(fileName)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()


def saveEtreeToFile(etree_element, fileName: Union[str, Path], pretty_print: bool = True):
    """
    Save an lxml ElementTree element to an XML file.

    Legacy compatibility function used by plugins to write XML output.
    Supports both lxml and xml.etree.ElementTree formats.

    Args:
        etree_element: lxml or ElementTree element to save
        fileName: Path to output file
        pretty_print: If True, format with indentation (default: True)

    Raises:
        IOError: If file cannot be written

    Example:
        >>> from lxml import etree
        >>> root = etree.Element('root')
        >>> saveEtreeToFile(root, '/path/to/output.xml')
    """
    file_path = Path(fileName)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Try lxml first (preferred for legacy compatibility)
    try:
        from lxml import etree as lxml_etree
        # Check if this is an lxml element
        if isinstance(etree_element, lxml_etree._Element):
            xml_bytes = lxml_etree.tostring(
                etree_element,
                pretty_print=pretty_print,
                xml_declaration=True,
                encoding='utf-8'
            )
            with open(file_path, 'wb') as f:
                f.write(xml_bytes)
            return
    except ImportError:
        pass

    # Fall back to standard library ElementTree
    import xml.etree.ElementTree as ET

    # If pretty_print is requested, indent the tree
    if pretty_print:
        try:
            ET.indent(etree_element)  # Python 3.9+
        except AttributeError:
            # For Python < 3.9, write without indenting
            pass

    tree = ET.ElementTree(etree_element)
    tree.write(str(file_path), encoding='utf-8', xml_declaration=True)


def backupFile(file_path: Union[str, Path], delete: bool = True) -> Optional[Path]:
    """
    Backup a file by renaming it with a numeric suffix.

    Args:
        file_path: Path to file to backup
        delete: If True, delete the original file after backup. If False, leave it in place.

    Returns:
        Path to backup file if created, None if original file doesn't exist

    Example:
        >>> backupFile('/path/to/file.txt', delete=False)
        Path('/path/to/file.txt.1')
    """
    import shutil

    file_path = Path(file_path)

    if not file_path.exists():
        return None

    # Find next available backup number
    counter = 1
    while True:
        backup_path = file_path.with_suffix(file_path.suffix + f'.{counter}')
        if not backup_path.exists():
            break
        counter += 1

    # Create backup
    if delete:
        # Move (rename) the file
        shutil.move(str(file_path), str(backup_path))
    else:
        # Copy the file
        shutil.copy2(str(file_path), str(backup_path))

    return backup_path
