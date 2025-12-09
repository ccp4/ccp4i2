"""
Debug why get_leaf_paths misses NCYCLES during traversal.
"""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "server"))

# Set CCP4I2_ROOT for plugin discovery
os.environ["CCP4I2_ROOT"] = str(project_root)

from ccp4i2.core.CCP4Modules import TASKMANAGER
from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core import CCP4Container as CCP4Container_module


def traverse_debug(node, path_parts, depth=0):
    """Debug version of traverse that prints what it's doing."""
    indent = "  " * depth
    leaf_paths = []

    print(f"{indent}Traversing: {node.objectName()} (type={type(node).__name__}, children={len(list(node.children()))})")

    for i, child in enumerate(node.children()):
        child_name = child.objectName() if hasattr(child, 'objectName') else str(child)
        child_type = type(child).__name__

        if isinstance(child, (CContainer, CCP4Container_module.CContainer)):
            print(f"{indent}  [{i}] {child_name} -> CONTAINER, recursing...")
            leaf_paths.extend(traverse_debug(child, path_parts + [child_name], depth + 1))
        else:
            print(f"{indent}  [{i}] {child_name} -> LEAF ({child_type})")
            qualifiers_dict = getattr(child, 'qualifiers', {})
            leaf_paths.append(
                {
                    "path": child.objectPath(),
                    "minimumPath": child.objectPath(),
                    "qualifiers": qualifiers_dict if isinstance(qualifiers_dict, dict) else {},
                    "className": type(child).__name__,
                    "object": child,
                }
            )

    return leaf_paths


def test_traverse_control_parameters():
    """Test traversing just controlParameters."""

    # Get plugin
    plugin_class = TASKMANAGER().get_plugin_class('prosmart_refmac')
    plugin = plugin_class(parent=None, workDirectory=None)

    ctrl = plugin.container.controlParameters

    print(f"=== Traversing controlParameters ===")
    print(f"controlParameters has {len(list(ctrl.children()))} children")
    print(f"NCYCLES is accessible: {hasattr(ctrl, 'NCYCLES')}")
    print(f"NCYCLES value: {ctrl.NCYCLES.value if hasattr(ctrl, 'NCYCLES') else 'N/A'}")

    print(f"\n=== Starting traverse ===")
    leaf_paths = traverse_debug(ctrl, [ctrl.objectName()], depth=0)

    print(f"\n=== Results ===")
    print(f"Found {len(leaf_paths)} leaves")

    ncycles_found = any('NCYCLES' in leaf['path'] for leaf in leaf_paths)
    print(f"NCYCLES found: {ncycles_found}")

    if not ncycles_found:
        print("\n=== Direct children inspection ===")
        for i, child in enumerate(ctrl.children()):
            child_name = child.objectName() if hasattr(child, 'objectName') else str(child)
            if 'NCYCLES' in child_name or i < 5:  # Show first 5 and any NCYCLES
                print(f"  [{i}] {child_name}: {type(child).__name__}")


if __name__ == "__main__":
    test_traverse_control_parameters()
