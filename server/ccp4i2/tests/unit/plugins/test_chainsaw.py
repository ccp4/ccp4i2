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
Tests for building and using chainsaw.
"""

import pytest
from ccp4i2 import I2_TOP
from ccp4i2.core.tasks import get_plugin_class


def test_chainsaw():
    task = get_plugin_class("chainsaw")()
    task.container.inputData.XYZIN = str(I2_TOP / "demo_data" / "mdm2" / "4hg7.pdb")
    assert task.container.inputData.XYZIN.__str__().endswith("4hg7.pdb")
