"""A CData object serialises its declared data-attributes and nothing else.

The fileContent property on CDataFile auto-creates a content object parented
to the file (so the hierarchy works), under the name 'fileContent'. It is a
runtime cache of what is inside the file --- not a data-attribute --- yet
dataOrder()'s remaining-children fallback used to sweep it into params.xml:
a gleaned CSeqDataFile carried its whole sequence text, a CAsuDataFile its
whole seqList, alongside the real attributes (baseName, relPath, dbFileId).

The fallback exists for CONTAINMENT --- CContainer children come from
def.xml, not class metadata, and must keep serialising. COMPOSITION classes
(declared content() fields) now serialise only what they declare.

Pure Python --- no CCP4 binaries.
"""
import pathlib
import xml.etree.ElementTree as ET

import pytest

from ccp4i2.core import CCP4ModelData

FASTA = (pathlib.Path(__file__).parents[3] / 'demo_data'
         / 'CDK1CyclinBCKS2' / 'P06493.fasta')


def _loaded_seq_file():
    f = CCP4ModelData.CSeqDataFile()
    f.setFullPath(str(FASTA))
    # Touch the property: auto-creates and loads the content object, which
    # becomes a hierarchy child named 'fileContent'
    content = f.fileContent
    assert content is not None
    return f


def test_loaded_file_content_is_a_child_but_not_serialised():
    f = _loaded_seq_file()
    child_names = {c.objectName() for c in f.children()
                   if hasattr(c, 'objectName')}
    assert 'fileContent' in child_names, \
        'precondition gone: fileContent no longer parents to the file, ' \
        'so this test guards nothing'

    elem = f.getEtree('SEQIN')
    tags = {e.tag for e in elem.iter()}
    assert 'fileContent' not in tags, \
        'the loaded content object leaked into the serialisation'
    assert 'sequence' not in tags, \
        'file content fields leaked into the serialisation'


def test_the_declared_attributes_still_serialise():
    f = _loaded_seq_file()
    elem = f.getEtree('SEQIN')
    tags = {e.tag for e in elem.iter()}
    assert {'baseName', 'relPath'} <= tags, \
        f'declared attributes went missing: {sorted(tags)}'
    base = elem.find('baseName')
    assert base is not None and base.text.endswith('P06493.fasta')


def test_data_order_lists_only_declared_names():
    f = _loaded_seq_file()
    order = f.dataOrder()
    assert 'fileContent' not in order
    assert 'baseName' in order


def test_container_children_still_serialise():
    """The fallback must survive for containment: a CContainer's children
    come from def.xml, not class metadata."""
    from ccp4i2.core.base_object.ccontainer import CContainer
    from ccp4i2.core.base_object.fundamental_types import CInt

    c = CContainer(name='parent')
    c.addContent(CInt, 'NCYCLES')
    c.NCYCLES.set(5)
    elem = c.getEtree('parent')
    assert elem.find('NCYCLES') is not None, \
        'container serialisation lost its def.xml-style children'
