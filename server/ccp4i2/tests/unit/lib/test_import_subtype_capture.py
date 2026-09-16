"""Import-by-use captures a file with the slot's declared subtype.

A raw map browsed or fetched into a half-map input has no intrinsic subtype, so
the importer must fall back to the slot's requiredSubType (5) rather than a
blanket 1 -- otherwise the file is stored as an ordinary map and never
recognised as a half map downstream. This pins the primary-subtype helper.
"""

import pytest

from ccp4i2.lib.utils.files.upload_param import _primary_required_subtype


@pytest.mark.parametrize("value,expected", [
    ([5, 1, 0], 5),          # def_xml_handler list form -> primary is the intent
    ([5], 5),
    ("5,1,0", 5),            # legacy comma-string form
    (5, 5),
    ("4", 4),                # mask slot
    ([1, 0], 1),             # primary 1 -> historical default, unchanged
    ([0], 1),                # 0 is "no specific type" -> 1
    (0, 1),
    (None, 1),               # no requiredSubType at all -> 1
    ("", 1),
    ([], 1),
    ("A,D", 1),              # unparseable -> 1, never raises
])
def test_primary_required_subtype(value, expected):
    assert _primary_required_subtype(value) == expected
