"""
Test get_leaf_paths to debug why NCYCLES isn't appearing in keywords.
"""

from ccp4i2.core.CCP4Modules import TASKMANAGER


# Import get_leaf_paths - need to avoid Django imports
def get_leaf_paths_simple(container):
    """Simplified version without Django dependencies."""
    from ccp4i2.core.CCP4Container import CContainer
    from ccp4i2.core import CCP4Container as CCP4Container_module

    def traverse(node, path_parts):
        leaf_paths = []
        for child in node.children():
            if isinstance(child, (CContainer, CCP4Container_module.CContainer)):
                leaf_paths.extend(traverse(child, path_parts + [child.objectName()]))
            else:
                # Get qualifiers dict from modern CData (attribute, not method)
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

    all_leafs = traverse(container, [container.objectName()])
    # Filter so only one item per unique path is returned
    unique = {}
    for leaf in all_leafs:
        if leaf["path"] not in unique:
            unique[leaf["path"]] = leaf
    return list(unique.values())


def test_prosmart_refmac_leaf_paths():
    """Test that get_leaf_paths finds NCYCLES."""

    # Get plugin class
    plugin_class = TASKMANAGER().get_plugin_class('prosmart_refmac')

    # Instantiate plugin
    plugin = plugin_class(parent=None, workDirectory=None)

    print(f"Plugin: {plugin}")
    print(f"Container: {plugin.container}")
    print(f"Container has controlParameters: {hasattr(plugin.container, 'controlParameters')}")

    # Get leaf paths
    leaf_paths = get_leaf_paths_simple(plugin.container)

    print(f"\n=== Found {len(leaf_paths)} leaf paths ===")

    # Look for NCYCLES
    ncycles_found = False
    for leaf in leaf_paths:
        if 'NCYCLES' in leaf['path']:
            ncycles_found = True
            print(f"\n✓ Found NCYCLES:")
            print(f"  path: {leaf['path']}")
            print(f"  className: {leaf['className']}")
            print(f"  qualifiers: {leaf.get('qualifiers', {})}")
            break

    if not ncycles_found:
        print("\n✗ NCYCLES not found in leaf paths!")
        print("\n=== All paths containing 'control' ===")
        for leaf in leaf_paths:
            if 'control' in leaf['path'].lower():
                print(f"  {leaf['path']}")

        print("\n=== First 30 paths ===")
        for leaf in leaf_paths[:30]:
            print(f"  {leaf['path']}")

        # Try to directly check controlParameters
        print("\n=== Checking plugin.container.controlParameters directly ===")
        if hasattr(plugin.container, 'controlParameters'):
            ctrl = plugin.container.controlParameters
            print(f"controlParameters: {ctrl}")
            print(f"controlParameters children count: {len(list(ctrl.children()))}")
            print(f"controlParameters has NCYCLES: {hasattr(ctrl, 'NCYCLES')}")

            if hasattr(ctrl, 'NCYCLES'):
                print(f"NCYCLES object: {ctrl.NCYCLES}")
                print(f"NCYCLES objectPath: {ctrl.NCYCLES.objectPath()}")
                print(f"NCYCLES type: {type(ctrl.NCYCLES)}")
                print(f"NCYCLES is in children: {ctrl.NCYCLES in list(ctrl.children())}")

    assert ncycles_found, "NCYCLES should be found in leaf paths"


if __name__ == "__main__":
    try:
        test_prosmart_refmac_leaf_paths()
        print("\n✅ Test passed!")
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
