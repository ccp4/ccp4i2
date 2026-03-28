#!/usr/bin/env python3
# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Test CAsuContentSeqList serialization using i2run-style construction.
This mimics how i2run creates items and uses CList.set() to copy them.
"""
import xml.etree.ElementTree as ET

from ccp4i2.core.CCP4ModelData import CAsuContentSeq, CAsuContentSeqList


def test_clist_i2run_style_construction():
    """Test that CAsuContentSeqList items created via setattr serialize correctly after set()."""
    # Create source list (like inputData.ASU_CONTENT in i2run)
    source_list = CAsuContentSeqList(name="ASU_CONTENT")

    # Create first item using setattr (like i2run does)
    item1 = CAsuContentSeq()
    setattr(item1, 'sequence', 'HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW')
    setattr(item1, 'nCopies', 1)
    setattr(item1, 'name', 'BETA')
    setattr(item1, 'polymerType', 'PROTEIN')
    setattr(item1, 'description', 'Beta-lactamase')

    source_list.append(item1)
    assert len(source_list) == 1

    # Create second item
    item2 = CAsuContentSeq()
    setattr(item2, 'name', 'BLIP')
    setattr(item2, 'sequence', 'AGVMTGAKFTQIQFGMTRQQVLDIAGAENCETGGSFGDSIHCRGHAAGDYYAYATFGFTSAAADAKVDSKSQEKLLAPSAPTLTLAKFNQVTVGMTRAQVLATVGQGSCTTWSEYYPAYPSTAGVTLSLSCFDVDGYSSTGFYRGSAHLWFTDGVLQGKRQWDLV')
    setattr(item2, 'polymerType', 'PROTEIN')
    setattr(item2, 'nCopies', 1)
    setattr(item2, 'description', 'Beta-lactamase inhibitory protein')

    source_list.append(item2)
    assert len(source_list) == 2

    # Now create target list and use .set() to copy (like ProvideAsuContents does)
    target_list = CAsuContentSeqList(name="seqList")
    target_list.set(source_list)
    assert len(target_list) == 2

    # Verify items in target list have content
    for i, item in enumerate(target_list):
        if hasattr(item, 'name'):
            name_val = item.name.value if hasattr(item.name, 'value') else item.name
            assert name_val is not None, f"Item {i} name should not be None"
        if hasattr(item, 'sequence'):
            seq_val = item.sequence.value if hasattr(item.sequence, 'value') else item.sequence
            assert seq_val is not None and len(seq_val) > 0, f"Item {i} sequence should not be empty"

    # Serialize target list to XML and verify all items have content
    etree = target_list.getEtree()
    items = list(etree)
    assert len(items) == 2, f"Expected 2 CAsuContentSeq elements in XML, got {len(items)}"

    for i, item_elem in enumerate(items):
        children = list(item_elem)
        assert len(children) > 0, f"Item {i} is EMPTY in XML - no child elements"
