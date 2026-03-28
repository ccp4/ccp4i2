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
"""Test CAsuDataFile.saveFile() after using CList.set()"""

import tempfile
from pathlib import Path

from ccp4i2.core.CCP4ModelData import CAsuContentSeq, CAsuContentSeqList, CAsuDataFile


def test_asu_savefile_with_clist_set():
    """Verify that saveFile() preserves name elements after CList.set()."""
    # Create source list with items (like inputData.ASU_CONTENT)
    source_list = CAsuContentSeqList(name="ASU_CONTENT")

    item1 = CAsuContentSeq()
    item1.sequence = 'HPETLVKVK'
    item1.nCopies = 1
    item1.name = 'BETA'
    item1.polymerType = 'PROTEIN'
    item1.description = 'Beta-lactamase'
    source_list.append(item1)

    item2 = CAsuContentSeq()
    item2.name = 'BLIP'
    item2.sequence = 'AGVMTGAKFT'
    item2.polymerType = 'PROTEIN'
    item2.nCopies = 1
    item2.description = 'BLIP protein'
    source_list.append(item2)

    assert len(source_list) == 2
    assert source_list[0].isSet('name')
    assert source_list[1].isSet('name')

    # Create ASU data file and copy items using set()
    asu_file = CAsuDataFile()

    with tempfile.TemporaryDirectory() as tmpdir:
        test_path = Path(tmpdir) / "test.asu.xml"
        asu_file.relPath = str(tmpdir)
        asu_file.baseName = "test.asu.xml"

        asu_file.fileContent.seqList.set(source_list)
        assert len(asu_file.fileContent.seqList) == 2

        # Verify items are populated
        for i, item in enumerate(asu_file.fileContent.seqList):
            assert item.isSet('name')
            assert item.isSet('sequence')

        # Save and read back
        asu_file.saveFile()
        content = test_path.read_text()

        assert '<name>BETA</name>' in content
        assert '<name>BLIP</name>' in content
