#!/usr/bin/env python3
"""Test script to verify fileContent property works"""
import os
import sys
from pathlib import Path

# Set CCP4I2_ROOT environment variable using relative path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
os.environ['CCP4I2_ROOT'] = str(PROJECT_ROOT)
sys.path.insert(0, str(PROJECT_ROOT))

from ccp4i2.core.CCP4ModelData import CPdbDataFile

print("[TEST] Creating CPdbDataFile...")
pdb_file = CPdbDataFile(name="XYZOUT")

print(f"\n[TEST] pdb_file: {pdb_file}")
print(f"[TEST] pdb_file.content: {pdb_file.content}")
print(f"[TEST] fileContentClassName qualifier: {pdb_file.get_qualifier('fileContentClassName')}")

print(f"\n[TEST] Accessing pdb_file.fileContent...")
file_content = pdb_file.fileContent
print(f"[TEST] pdb_file.fileContent: {file_content}")

if file_content is None:
    print("\n[TEST] FAILED! fileContent is None")
else:
    print(f"\n[TEST] SUCCESS! fileContent is {type(file_content).__name__}")
    print(f"[TEST] Has loadFile method: {hasattr(file_content, 'loadFile')}")
