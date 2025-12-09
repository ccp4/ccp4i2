"""
Tests for building and using chainsaw.
"""

import pytest
import os
from ccp4i2.core.CCP4TaskManager import TASKMANAGER


@pytest.mark.skipif(
    'CCP4I2_ROOT' not in os.environ,
    reason="CCP4I2_ROOT environment variable not set"
)
def test_chainsaw():
    task = TASKMANAGER().get_plugin_class("chainsaw")()
    task.container.inputData.XYZIN = os.path.join(os.environ["CCP4I2_ROOT"], "demo_data", "mdm2", "4hg7.pdb")
    assert task.container.inputData.XYZIN.__str__().endswith("4hg7.pdb")
