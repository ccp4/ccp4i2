"""Standalone test for JSON encoder to debug inputData.XYZIN issue."""
import json
import os
import sys
from pathlib import Path

# Set up environment
os.environ.setdefault("CCP4I2_ROOT", str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "server"))
sys.path.insert(0, str(Path(__file__).parent.parent / "pipelines/prosmart_refmac/script"))

from core import CCP4TaskManager
from core import CCP4Container
from core.base_object import CData
from core.base_object.fundamental_types import CString
from server.ccp4x.lib.utils.containers.json_encoder import CCP4i2JsonEncoder
from prosmart_refmac import prosmart_refmac


def test_cdatafile_children_serialization():
    """Test that CDataFile children are properly serialized."""
    # Create plugin instance which loads container from .def.xml
    print("\n=== Creating prosmart_refmac plugin ===")
    plugin = prosmart_refmac()
    container = plugin.container

    # Compare inputData.XYZIN (CPdbDataFile) vs inputData.F_SIGF (CObsDataFile)
    input_data = container.inputData
    output_data = container.outputData

    print("\n=== Comparing XYZIN (CPdbDataFile) vs F_SIGF (CObsDataFile) ===")

    xyzin = input_data.XYZIN
    f_sigf = input_data.F_SIGF

    print(f"\nXYZIN type: {type(xyzin).__name__}")
    print(f"F_SIGF type: {type(f_sigf).__name__}")

    print(f"\nXYZIN children ({len(xyzin.children())}):")
    for child in xyzin.children():
        name = child.objectName() if hasattr(child, 'objectName') else str(child)
        print(f"  - {name}: {type(child).__name__}")

    print(f"\nF_SIGF children ({len(f_sigf.children())}):")
    for child in f_sigf.children():
        name = child.objectName() if hasattr(child, 'objectName') else str(child)
        print(f"  - {name}: {type(child).__name__}")

    # Set values on both
    print("\n=== Setting values on both ===")
    xyzin.baseName.set("test_model.pdb")
    xyzin.project.set("test-project-uuid")
    xyzin.relPath.set("CCP4_JOBS/job_1")

    f_sigf.baseName.set("test_data.mtz")
    f_sigf.project.set("test-project-uuid")
    f_sigf.relPath.set("CCP4_JOBS/job_1")

    # Serialize both
    print("\n=== Serializing XYZIN ===")
    xyzin_json = json.loads(json.dumps(xyzin, cls=CCP4i2JsonEncoder))
    print(f"XYZIN _value keys: {list(xyzin_json.get('_value', {}).keys())}")
    print(f"XYZIN baseName._value: {xyzin_json.get('_value', {}).get('baseName', {}).get('_value')}")

    print("\n=== Serializing F_SIGF ===")
    f_sigf_json = json.loads(json.dumps(f_sigf, cls=CCP4i2JsonEncoder))
    print(f"F_SIGF _value keys: {list(f_sigf_json.get('_value', {}).keys())}")
    print(f"F_SIGF baseName._value: {f_sigf_json.get('_value', {}).get('baseName', {}).get('_value')}")

    # Also check outputData.XYZOUT (also CPdbDataFile)
    print("\n=== Checking outputData.XYZOUT ===")
    xyzout = output_data.XYZOUT
    print(f"XYZOUT type: {type(xyzout).__name__}")
    print(f"XYZOUT children ({len(xyzout.children())}):")
    for child in xyzout.children():
        name = child.objectName() if hasattr(child, 'objectName') else str(child)
        print(f"  - {name}: {type(child).__name__}")

    xyzout.baseName.set("output_model.pdb")
    xyzout_json = json.loads(json.dumps(xyzout, cls=CCP4i2JsonEncoder))
    print(f"XYZOUT _value keys: {list(xyzout_json.get('_value', {}).keys())}")
    print(f"XYZOUT baseName._value: {xyzout_json.get('_value', {}).get('baseName', {}).get('_value')}")

    print(f"\n=== XYZIN Type Info ===")
    print(f"XYZIN type: {type(xyzin)}")
    print(f"XYZIN class name: {xyzin.__class__.__name__}")
    print(f"XYZIN MRO: {[c.__name__ for c in type(xyzin).__mro__[:5]]}")

    # Check children
    children = xyzin.children()
    print(f"\n=== XYZIN Children ({len(children)}) ===")
    for child in children:
        child_name = child.objectName() if hasattr(child, 'objectName') else str(child)
        child_type = type(child).__name__
        child_value = getattr(child, '_value', 'NO_VALUE_ATTR')
        print(f"  - {child_name}: {child_type} = {repr(child_value)[:50]}")

    # Check if there are attributes that aren't children
    print(f"\n=== XYZIN __dict__ keys ===")
    for key in sorted(xyzin.__dict__.keys()):
        if not key.startswith('_'):
            val = xyzin.__dict__[key]
            print(f"  - {key}: {type(val).__name__}")

    # Check specific attributes via direct access
    print(f"\n=== Direct Attribute Access ===")
    for attr_name in ['baseName', 'project', 'relPath', 'annotation', 'dbFileId', 'subType', 'contentFlag', 'selection']:
        if hasattr(xyzin, attr_name):
            attr = getattr(xyzin, attr_name)
            attr_value = getattr(attr, '_value', 'NO_VALUE_ATTR') if attr else None
            is_child = any(c.objectName() == attr_name for c in children if hasattr(c, 'objectName'))
            print(f"  {attr_name}: {type(attr).__name__ if attr else 'None'}, _value={repr(attr_value)[:30]}, in_children={is_child}")
        else:
            print(f"  {attr_name}: NOT PRESENT")

    # Set some values to test
    print(f"\n=== Setting Values ===")
    xyzin.baseName.set("test_model.pdb")
    xyzin.project.set("test-project-uuid")
    xyzin.relPath.set("CCP4_JOBS/job_1")

    print(f"baseName after set: {repr(xyzin.baseName._value)}")
    print(f"project after set: {repr(xyzin.project._value)}")
    print(f"relPath after set: {repr(xyzin.relPath._value)}")

    # Check children again after setting values
    children_after = xyzin.children()
    print(f"\n=== XYZIN Children After Setting Values ({len(children_after)}) ===")
    for child in children_after:
        child_name = child.objectName() if hasattr(child, 'objectName') else str(child)
        child_type = type(child).__name__
        child_value = getattr(child, '_value', 'NO_VALUE_ATTR')
        print(f"  - {child_name}: {child_type} = {repr(child_value)[:50]}")

    # Serialize and check
    print(f"\n=== JSON Serialization ===")
    json_str = json.dumps(xyzin, cls=CCP4i2JsonEncoder, indent=2)
    json_obj = json.loads(json_str)

    print(f"_class: {json_obj.get('_class')}")
    print(f"_baseClass: {json_obj.get('_baseClass')}")
    print(f"_CONTENTS_ORDER: {json_obj.get('_CONTENTS_ORDER')}")

    xyzin_value = json_obj.get('_value', {})
    print(f"\n_value keys: {list(xyzin_value.keys())}")

    for key in ['baseName', 'project', 'relPath']:
        if key in xyzin_value:
            print(f"  {key}._value = {repr(xyzin_value[key].get('_value'))}")
        else:
            print(f"  {key}: MISSING from _value!")

    # Assert
    assert 'baseName' in xyzin_value, "baseName should be in serialized _value"
    assert xyzin_value['baseName'].get('_value') == 'test_model.pdb', \
        f"Expected 'test_model.pdb', got {repr(xyzin_value['baseName'].get('_value'))}"

    print("\n=== TEST PASSED ===")


if __name__ == "__main__":
    test_cdatafile_children_serialization()
