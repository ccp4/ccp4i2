"""Import a sequence file (FASTA / PIR / plain) into the project."""
from pathlib import Path

from ccp4i2.wrappers.import_common import CImportFileBase


class ImportSequence(CImportFileBase):

    TASKNAME = 'ImportSequence'
    INPUT_PARAM = 'SEQIN'
    OUTPUT_PARAM = 'SEQOUT'
    ANNOTATION_PREFIX = 'Imported sequence'

    _FORMAT_BY_EXT = {
        '.fasta': 'fasta', '.fa': 'fasta', '.seq': 'fasta', '.pir': 'pir',
    }

    def validate_source(self, src_path):
        suffix = Path(src_path).suffix.lower()
        fmt = self._FORMAT_BY_EXT.get(suffix, 'fasta')
        try:
            from Bio import SeqIO
            records = list(SeqIO.parse(src_path, fmt))
            if not records:
                # Fall back to a plain residue-string read
                with open(src_path) as handle:
                    body = handle.read().strip()
                if not body:
                    return 'Sequence file is empty.'
                return None
        except Exception as e:
            return 'Not a readable sequence file: ' + str(e)
        return None
