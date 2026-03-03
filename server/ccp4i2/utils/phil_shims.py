"""
PHIL Shims - Convert rich CCP4i2 data types to PHIL-compatible values.

Shims bridge the gap between CCP4i2's rich file types (CMtzDataFile,
CPdbDataFile, CAsuDataFile, CDictDataFile) and the simpler types that
PHIL understands (file paths as strings). They run at task execution
time, converting CCP4i2 input data into (phil_path, value) pairs that
get merged into the working_phil.

Usage:
    shims = [
        MtzFileShim("HKLIN", "picard.hklin"),
        PdbFileShim("XYZIN", ["picard.xyzin"]),
        AsuContentShim("ASUIN", "picard.seqin"),
        DictFileShim("DICT", "phasertng.cluster_compound.filename"),
    ]

    for shim in shims:
        params.extend(shim.convert(container, work_directory))
"""

import os
import logging

logger = logging.getLogger(__name__)


class PhilShim:
    """Base class for converting CCP4i2 data objects to PHIL parameter values."""

    def convert(self, container, work_directory):
        """Convert CCP4i2 data to PHIL parameter name=value pairs.

        Args:
            container: The plugin's container (with inputData, outputData, etc.)
            work_directory: The job's working directory (for writing intermediate files)

        Returns:
            list of (phil_path, value) tuples
        """
        raise NotImplementedError


class MtzFileShim(PhilShim):
    """Convert CMtzDataFile to a PHIL hklin file path parameter."""

    def __init__(self, input_name, phil_hklin_path):
        """
        Args:
            input_name: Name of the CMtzDataFile in container.inputData (e.g., "HKLIN")
            phil_hklin_path: PHIL dotted path for the file (e.g., "picard.hklin")
        """
        self.input_name = input_name
        self.phil_hklin_path = phil_hklin_path

    def convert(self, container, work_directory):
        result = []
        try:
            mtz_obj = getattr(container.inputData, self.input_name)
        except AttributeError:
            return result

        if mtz_obj is not None and mtz_obj.isSet():
            file_path = mtz_obj.getFullPath()
            if file_path:
                result.append((self.phil_hklin_path, file_path))
        return result


class PdbFileShim(PhilShim):
    """Convert CPdbDataFile to one or more PHIL file path parameters."""

    def __init__(self, input_name, phil_paths):
        """
        Args:
            input_name: Name of the CPdbDataFile in container.inputData (e.g., "XYZIN")
            phil_paths: PHIL path or list of paths to set (e.g., ["picard.xyzin"])
        """
        self.input_name = input_name
        if isinstance(phil_paths, str):
            phil_paths = [phil_paths]
        self.phil_paths = phil_paths

    def convert(self, container, work_directory):
        result = []
        try:
            pdb_obj = getattr(container.inputData, self.input_name)
        except AttributeError:
            return result

        if pdb_obj is not None and pdb_obj.isSet():
            file_path = pdb_obj.getFullPath()
            if file_path:
                for phil_path in self.phil_paths:
                    result.append((phil_path, file_path))
        return result


class PdbFileListShim(PhilShim):
    """Convert a CList of CPdbDataFile to repeated PHIL file path parameters.

    For PHIL parameters defined with `.multiple = True`, each list item
    produces a separate (phil_path, value) tuple with the same path name.
    PHIL's fetch() merges these into the repeated parameter correctly.
    """

    def __init__(self, input_name, phil_path):
        """
        Args:
            input_name: Name of the CList in container.inputData (e.g., "XYZIN")
            phil_path: PHIL dotted path for each file (e.g., "picard.xyzin")
        """
        self.input_name = input_name
        self.phil_path = phil_path

    def convert(self, container, work_directory):
        result = []
        try:
            list_obj = getattr(container.inputData, self.input_name)
        except AttributeError:
            return result

        for item in list_obj:
            if item is not None and item.isSet():
                file_path = item.getFullPath()
                if file_path:
                    result.append((self.phil_path, file_path))
        return result


class AsuContentShim(PhilShim):
    """Convert CAsuDataFile to FASTA file(s) and PHIL sequence path parameters.

    This is the most complex shim. CAsuDataFile wraps CAsuContent, which holds
    a seqList of CAsuContentSeq objects. Each sequence has:
        - sequence: the amino acid / nucleotide string
        - nCopies: number of copies in the ASU
        - polymerType: "PROTEIN", "RNA", "DNA"
        - name: identifier
        - description: human-readable description

    The shim writes these sequences to FASTA format and returns the path.
    """

    def __init__(self, input_name, phil_seq_path):
        """
        Args:
            input_name: Name of the CAsuDataFile in container.inputData (e.g., "ASUIN")
            phil_seq_path: PHIL dotted path for the sequence file (e.g., "picard.seqin")
        """
        self.input_name = input_name
        self.phil_seq_path = phil_seq_path

    def convert(self, container, work_directory):
        result = []
        try:
            asu_obj = getattr(container.inputData, self.input_name)
        except AttributeError:
            return result

        if asu_obj is None or not asu_obj.isSet():
            return result

        # Load the ASU file content
        try:
            asu_obj.loadFile()
        except Exception as e:
            logger.warning("Could not load ASU file for %s: %s", self.input_name, e)
            return result

        # Get the file content (CAsuContent)
        content = getattr(asu_obj, 'fileContent', None)
        if content is None:
            # Try accessing via fileContent attribute (some implementations)
            content = asu_obj
            if not hasattr(content, 'seqList'):
                return result

        seq_list = getattr(content, 'seqList', None)
        if seq_list is None:
            return result

        # Write sequences to FASTA file
        fasta_path = os.path.join(work_directory, "sequences.fasta")
        sequences_written = 0

        try:
            with open(fasta_path, "w") as f:
                for seq_obj in seq_list:
                    # Extract sequence string
                    sequence = seq_obj.sequence
                    if hasattr(sequence, 'value'):
                        sequence = str(sequence.value)
                    else:
                        sequence = str(sequence)

                    # Clean whitespace
                    sequence = sequence.replace(' ', '').replace('\n', '').replace('\t', '')

                    if not sequence:
                        continue

                    # Extract name
                    name = getattr(seq_obj, 'name', None)
                    if name is not None and hasattr(name, 'value'):
                        name = str(name.value)
                    elif name is not None:
                        name = str(name)
                    else:
                        name = f"sequence_{sequences_written + 1}"

                    # Write FASTA entry
                    f.write(f">{name}\n")
                    # Wrap at 80 characters
                    for i in range(0, len(sequence), 80):
                        f.write(sequence[i:i + 80] + "\n")

                    sequences_written += 1

        except Exception as e:
            logger.error("Error writing FASTA file: %s", e)
            return result

        if sequences_written > 0:
            result.append((self.phil_seq_path, fasta_path))

        return result


class DictFileShim(PhilShim):
    """Convert CDictDataFile to a PHIL dictionary/compound file path parameter."""

    def __init__(self, input_name, phil_path):
        """
        Args:
            input_name: Name of the CDictDataFile in container.inputData (e.g., "DICT")
            phil_path: PHIL dotted path for the file
        """
        self.input_name = input_name
        self.phil_path = phil_path

    def convert(self, container, work_directory):
        result = []
        try:
            dict_obj = getattr(container.inputData, self.input_name)
        except AttributeError:
            return result

        if dict_obj is not None and dict_obj.isSet():
            file_path = dict_obj.getFullPath()
            if file_path:
                result.append((self.phil_path, file_path))
        return result


class DagFileShim(PhilShim):
    """Convert CPhaserTngDagFile to a PHIL DAG file path parameter."""

    def __init__(self, input_name, phil_path):
        """
        Args:
            input_name: Name of the CPhaserTngDagFile in container.inputData (e.g., "DAGIN")
            phil_path: PHIL dotted path for the file (e.g., "phasertng.put_solution.dag.filename")
        """
        self.input_name = input_name
        self.phil_path = phil_path

    def convert(self, container, work_directory):
        result = []
        try:
            dag_obj = getattr(container.inputData, self.input_name)
        except AttributeError:
            return result

        if dag_obj is not None and dag_obj.isSet():
            file_path = dag_obj.getFullPath()
            if file_path:
                result.append((self.phil_path, file_path))
        return result
