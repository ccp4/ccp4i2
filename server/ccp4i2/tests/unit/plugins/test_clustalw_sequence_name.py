"""ClustalW names a PDB sequence by its entry and chain, not "CHAIN".

The report names each sequence by the third |-separated field of its
header: right for UniProt, but every RCSB header reads
">4HG7:A|PDBID|CHAIN|SEQUENCE", so every PDB sequence was called "CHAIN".
"""
import pytest

from ccp4i2.wrappers.clustalw.script.clustalw import sequence_name


@pytest.mark.parametrize("header, name", [
    ("4HG7:A|PDBID|CHAIN|SEQUENCE", "4HG7_A"),
    ("4HG7_1|Chain A|E3 ubiquitin-protein ligase Mdm2|Homo sapiens", "4HG7_1"),
    ("sp|Q01094|E2F1_HUMAN Transcription factor E2F1", "sp|Q01094|E2F1_HUMAN"),
    ("MDMX 3dab chain A", "MDMX"),
    ("", ""),
])
def test_sequence_name(header, name):
    assert sequence_name(header) == name
