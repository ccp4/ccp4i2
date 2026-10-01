"""The sequence selection on an AU contents file reaches the program.

MolRep failed on a two-sequence AU file with one sequence selected
(issue #669): the selection CDict kept its entries where no generic CData
path looked, so it was never marked set, was written to params.xml as an
empty element, and was dropped by copyData into the pipeline's sub-job.
"""

import xml.etree.ElementTree as ET

import pytest

from ccp4i2.core.base_object.ccontainer import CContainer
from ccp4i2.core.CCP4Data import CDict
from ccp4i2.core.CCP4ModelData import CAsuContentSeq, CAsuDataFile


def _seq(name, sequence):
    seq = CAsuContentSeq()
    seq.name = name
    seq.sequence = sequence
    seq.polymerType = 'PROTEIN'
    seq.nCopies = 1
    return seq


@pytest.fixture
def two_seq_asu(tmp_path):
    """A saved AU contents file holding two sequences."""
    asu = CAsuDataFile()
    asu.relPath = str(tmp_path)
    asu.baseName = 'two.asu.xml'
    asu.fileContent.seqList.append(_seq('LMX', 'MKTAYIAKQRQISFVKSHFSRQ'))
    asu.fileContent.seqList.append(_seq('UBQ', 'MQIFVKTLTGKTITLEVEPSDT'))
    asu.saveFile()
    return tmp_path / 'two.asu.xml'


def _asu(path, mode, name='ASUIN', parent=None):
    asu = CAsuDataFile(name=name, parent=parent)
    asu.set_qualifier('selectionMode', mode)
    asu.setFullPath(str(path))
    return asu


def _fasta_names(asu, tmp_path):
    out = tmp_path / 'out.fasta'
    asu.writeFasta(str(out))
    return [line[1:].strip() for line in out.read_text().splitlines()
            if line.startswith('>')]


def _codes(report):
    return [e['code'] for e in report._errors]


# --- CDict -----------------------------------------------------------------

def test_cdict_is_set_once_it_has_entries():
    d = CDict(name='selection')
    assert not d.isSet()
    d.update({'LMX': True, 'UBQ': False})  # what set_parameter does
    assert d.isSet()
    d.unSet()
    assert not d.isSet() and len(d) == 0


def test_cdict_xml_round_trip_in_qt_layout():
    d = CDict(name='selection')
    d.set({'LMX': True, 'UBQ': False})
    elem = d.getEtree()
    items = [(i.findtext('key'), i.findtext('value')) for i in elem.findall('item')]
    assert items == [('LMX', 'True'), ('UBQ', 'False')]

    back = CDict(name='selection')
    back.setEtree(ET.fromstring(ET.tostring(elem)))
    assert dict(back.items()) == {'LMX': True, 'UBQ': False}


# --- writeFasta and copyData ---------------------------------------------

def test_write_fasta_honours_the_selection(two_seq_asu, tmp_path):
    asu = _asu(two_seq_asu, 1)
    asu.selection.update({'LMX': False, 'UBQ': True})
    assert _fasta_names(asu, tmp_path) == ['UBQ']


def test_write_fasta_ignores_a_selection_in_mode_0(two_seq_asu, tmp_path):
    asu = _asu(two_seq_asu, '0')  # def.xml qualifiers arrive as strings
    asu.selection.update({'LMX': False})
    assert _fasta_names(asu, tmp_path) == ['LMX', 'UBQ']


def test_copydata_carries_the_selection_to_a_subjob(two_seq_asu, tmp_path):
    pipe = CContainer(name='inputData')
    pipe.ASUIN = _asu(two_seq_asu, 1, parent=pipe)
    pipe._data_order.append('ASUIN')
    pipe.ASUIN.selection.update({'LMX': False, 'UBQ': True})

    sub = CContainer(name='inputData')
    sub.ASUIN = _asu(two_seq_asu, 1, parent=sub)
    sub._data_order.append('ASUIN')
    sub.ASUIN.selection.unSet()
    sub.copyData(pipe)

    assert dict(sub.ASUIN.selection.items()) == {'LMX': False, 'UBQ': True}
    assert _fasta_names(sub.ASUIN, tmp_path) == ['UBQ']


def test_selection_survives_the_params_xml_round_trip(two_seq_asu):
    asu = _asu(two_seq_asu, 1)
    asu.selection.update({'LMX': False, 'UBQ': True})
    elem = asu.getEtree()

    back = _asu(two_seq_asu, 1)
    back.setEtree(ET.fromstring(ET.tostring(elem)))
    assert dict(back.selection.items()) == {'LMX': False, 'UBQ': True}


def test_params_handler_reads_the_selection_back(two_seq_asu):
    """The job runner reloads params.xml through ParamsXmlHandler, whose
    attribute walk skipped a dict's <item> entries: the next save then
    wrote the selection away."""
    from ccp4i2.core.task_manager.params_xml_handler import ParamsXmlHandler

    asu = _asu(two_seq_asu, 1)
    asu.selection.update({'LMX': False, 'UBQ': True})
    elem = ET.fromstring(ET.tostring(asu.getEtree('ASUIN')))

    back = _asu(two_seq_asu, 1)
    back.selection.unSet()
    ParamsXmlHandler()._import_structured_data(elem, back)
    assert dict(back.selection.items()) == {'LMX': False, 'UBQ': True}


# --- validity: the selectionMode rule ------------------------------------

@pytest.mark.parametrize('mode,selection,codes', [
    (1, None, [110]),                                # nothing chosen: both count
    (1, {'LMX': True, 'UBQ': True}, [110]),
    (1, {'LMX': False, 'UBQ': True}, []),
    ('1', {'LMX': False, 'UBQ': False}, [110]),
    (2, None, []),
    (2, {'LMX': False, 'UBQ': False}, [111]),
    (0, {'LMX': False, 'UBQ': False}, []),
])
def test_validity_follows_selection_mode(two_seq_asu, mode, selection, codes):
    asu = _asu(two_seq_asu, mode)
    if selection:
        asu.selection.update(selection)
    assert _codes(asu.validity()) == codes


def test_one_sequence_needs_no_choice(tmp_path):
    asu = CAsuDataFile()
    asu.relPath = str(tmp_path)
    asu.baseName = 'one.asu.xml'
    asu.fileContent.seqList.append(_seq('LMX', 'MKTAYIAKQRQ'))
    asu.saveFile()
    assert _codes(_asu(tmp_path / 'one.asu.xml', 1).validity()) == []
