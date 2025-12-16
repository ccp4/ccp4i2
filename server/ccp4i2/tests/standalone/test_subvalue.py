#!/usr/bin/env python
"""Test script to verify subValue argument handling."""

import os
import sys
from pathlib import Path

# Set up environment using relative path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
os.environ['CCP4I2_ROOT'] = str(PROJECT_ROOT)
os.environ['DJANGO_SETTINGS_MODULE'] = 'ccp4i2.config.settings'

# Add server to path
sys.path.insert(0, str(PROJECT_ROOT / 'server'))

# Initialize Django
import django
django.setup()

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from ccp4i2.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase

# Test 1: fullPath subValue
print("=" * 60)
print("Test 1: Testing --XYZIN fullPath=/path/to/file.pdb")
print("=" * 60)

args1 = [
    "prosmart_refmac",
    "--XYZIN",
    "fullPath=/tmp/test.pdb",
]

try:
    runner1 = CCP4i2RunnerBase(the_args=args1)
    plugin1 = runner1.getPlugin()

    # Check if XYZIN.fullPath was set
    from ccp4i2.lib.utils.containers.find_objects import find_object_by_path
    xyzin = find_object_by_path(plugin1.container, "inputData.XYZIN")
    if xyzin:
        print(f"✓ XYZIN object found")
        print(f"  fullPath: {xyzin.fullPath}")
        print(f"  fullPath is set: {xyzin.isSet('fullPath')}")
    else:
        print("✗ XYZIN object not found")
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 2: composite path with /
print("\n" + "=" * 60)
print("Test 2: Testing --XYZIN selection/text=(NUT)")
print("=" * 60)

args2 = [
    "prosmart_refmac",
    "--XYZIN",
    "/tmp/test.pdb",
    "selection/text=(NUT)",
]

try:
    runner2 = CCP4i2RunnerBase(the_args=args2)
    plugin2 = runner2.getPlugin()

    # Check if XYZIN.selection.text was set
    from ccp4i2.lib.utils.containers.find_objects import find_object_by_path
    xyzin = find_object_by_path(plugin2.container, "inputData.XYZIN")
    if xyzin:
        print(f"✓ XYZIN object found")
        selection = find_object_by_path(xyzin, "selection")
        if selection:
            print("✓ selection object found")
            text_obj = getattr(selection, 'text', None)
            if text_obj:
                print(f"  selection.text.value: {text_obj.value if hasattr(text_obj, 'value') else text_obj}")
                print(f"  selection.text is set: {text_obj.isSet('value') if hasattr(text_obj, 'isSet') else 'N/A'}")
        else:
            print("✗ selection object not found")
    else:
        print("✗ XYZIN object not found")
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Combined - fullPath and selection
print("\n" + "=" * 60)
print("Test 3: Testing --XYZIN fullPath=/path selection/text=(NUT)")
print("=" * 60)

args3 = [
    "prosmart_refmac",
    "--XYZIN",
    "fullPath=/tmp/test.pdb",
    "selection/text=(NUT)",
]

try:
    runner3 = CCP4i2RunnerBase(the_args=args3)
    plugin3 = runner3.getPlugin()

    # Check if both were set
    from ccp4i2.lib.utils.containers.find_objects import find_object_by_path
    xyzin = find_object_by_path(plugin3.container, "inputData.XYZIN")
    if xyzin:
        print(f"✓ XYZIN object found")
        print(f"  fullPath: {xyzin.fullPath}")
        print(f"  fullPath is set: {xyzin.isSet('fullPath')}")

        selection = find_object_by_path(xyzin, "selection")
        if selection:
            print("✓ selection object found")
            text_obj = getattr(selection, 'text', None)
            if text_obj:
                print(f"  selection.text.value: {text_obj.value if hasattr(text_obj, 'value') else text_obj}")
                print(f"  selection.text is set: {text_obj.isSet('value') if hasattr(text_obj, 'isSet') else 'N/A'}")
        else:
            print("✗ selection object not found")
    else:
        print("✗ XYZIN object not found")
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 60)
print("All tests completed")
print("=" * 60)
