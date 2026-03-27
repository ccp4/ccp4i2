"""
Tests for building and using chainsaw.
"""

import pytest
from ccp4i2 import I2_TOP
from ccp4i2.core.tasks import get_plugin_class


def test_chainsaw():
    task = get_plugin_class("chainsaw")()
    task.container.inputData.XYZIN = str(I2_TOP / "demo_data" / "mdm2" / "4hg7.pdb")
    assert task.container.inputData.XYZIN.__str__().endswith("4hg7.pdb")
