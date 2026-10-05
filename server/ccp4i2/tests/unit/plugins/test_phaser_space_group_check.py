"""Phaser's space-group change, told apart from a renaming by number.

A held-out agent run saw P 2 21 21 become P 21 21 21 after MR, read it as a
setting difference, and was unsure which data to carry on with. They are
space groups 18 and 19.
"""
import pytest

pytest.importorskip("gemmi")

from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil import (  # noqa: E402
    space_group_check,
)


@pytest.mark.parametrize("given, solved, change", [
    ("P 2 21 21", "P 21 21 21", "changed"),   # 18 -> 19: a screw axis gained
    ("P 61 2 2", "P 65 2 2", "changed"),      # the enantiomorph Phaser chose
    ("P 1 21 1", "P 21", "setting"),          # one group, two names
    ("P 21 21 2", "P 2 21 21", "setting"),    # No. 18 in another setting
    ("P 21 21 21", "P 21 21 21", "none"),
])
def test_change_is_decided_by_space_group_number(given, solved, change):
    node = space_group_check(given, solved)
    assert node.get("change") == change
    assert node.get("given") == given and node.get("solved") == solved
