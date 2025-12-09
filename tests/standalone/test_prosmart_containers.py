#!/usr/bin/env python3
"""
Test script to verify that prosmart_refmac loads all containers from .def.xml
"""
import os
import sys
from pathlib import Path

# Set up environment using relative path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
os.environ['CCP4I2_ROOT'] = str(PROJECT_ROOT)

# Add project root to path
sys.path.insert(0, str(PROJECT_ROOT))

# Import the plugin class directly
sys.path.insert(0, str(PROJECT_ROOT / 'pipelines/prosmart_refmac/script'))
from prosmart_refmac import prosmart_refmac

print("[TEST] Creating prosmart_refmac plugin instance...")
plugin = prosmart_refmac()

print(f"\n[TEST] Plugin created: {plugin}")
print(f"[TEST] Plugin TASKNAME: {plugin.TASKNAME}")
print(f"[TEST] Plugin container: {plugin.container}")

print(f"\n[TEST] Checking for standard containers:")
print(f"  - inputData: {hasattr(plugin.container, 'inputData')}")
print(f"  - outputData: {hasattr(plugin.container, 'outputData')}")
print(f"  - controlParameters: {hasattr(plugin.container, 'controlParameters')}")
print(f"  - guiAdmin: {hasattr(plugin.container, 'guiAdmin')}")

print(f"\n[TEST] Checking for prosmart_refmac-specific containers:")
print(f"  - prosmartProtein: {hasattr(plugin.container, 'prosmartProtein')}")
print(f"  - prosmartNucleicAcid: {hasattr(plugin.container, 'prosmartNucleicAcid')}")
print(f"  - platonyzer: {hasattr(plugin.container, 'platonyzer')}")

print(f"\n[TEST] All children of plugin.container:")
for child in plugin.container.children():
    print(f"  - {child.name} ({type(child).__name__})")

if hasattr(plugin.container, 'prosmartProtein'):
    print(f"\n[TEST] SUCCESS! prosmartProtein container exists")
    print(f"[TEST] Checking prosmartProtein.TOGGLE...")
    if hasattr(plugin.container.prosmartProtein, 'TOGGLE'):
        print(f"  - TOGGLE exists: {plugin.container.prosmartProtein.TOGGLE}")
        print(f"  - TOGGLE value: {plugin.container.prosmartProtein.TOGGLE.value}")
    else:
        print(f"  - ERROR: TOGGLE not found in prosmartProtein!")
        print(f"  - prosmartProtein children:")
        for child in plugin.container.prosmartProtein.children():
            print(f"    - {child.name} ({type(child).__name__})")
else:
    print(f"\n[TEST] FAILED! prosmartProtein container does not exist")
    print(f"\n[TEST] Debugging: Let's check what the parser returned...")

    # Re-parse the .def.xml file to see what it actually contains
    from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
    parser = DefXmlParser()
    parsed = parser.parse_def_xml(str(PROJECT_ROOT / 'pipelines/prosmart_refmac/script/prosmart_refmac.def.xml'))
    print(f"\n[TEST] Parsed container has {len(parsed.children())} children:")
    for child in parsed.children():
        print(f"  - {child.name} ({type(child).__name__})")

    sys.exit(1)
