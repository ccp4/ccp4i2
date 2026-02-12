"""
Test that prosmart_refmac plugin keywords include inherited NCYCLES parameter.

This test verifies that when prosmart_refmac is instantiated as a plugin,
its container includes NCYCLES and other parameters inherited from refmac.
"""

from ccp4i2.core.task_manager.plugin_registry import get_plugin_class


def test_prosmart_refmac_plugin_has_ncycles():
    """Test that prosmart_refmac plugin instance has NCYCLES in container."""

    # Get plugin class
    plugin_class = get_plugin_class('prosmart_refmac')
    print(f"Plugin class: {plugin_class}")

    # Instantiate plugin
    plugin = plugin_class(parent=None, workDirectory=None)
    print(f"Plugin instance: {plugin}")
    print(f"Plugin container: {plugin.container}")

    # Check that container exists
    assert hasattr(plugin, 'container'), "Plugin should have container"

    # Check that controlParameters exists
    assert hasattr(plugin.container, 'controlParameters'), \
        "Plugin container should have controlParameters"

    # List all attributes in controlParameters
    print("\n=== controlParameters attributes ===")
    for attr_name in dir(plugin.container.controlParameters):
        if not attr_name.startswith('_'):
            print(f"  {attr_name}")

    # Check that NCYCLES exists
    assert hasattr(plugin.container.controlParameters, 'NCYCLES'), \
        "Plugin container.controlParameters should have NCYCLES (inherited from refmac)"

    print(f"\n✓ prosmart_refmac plugin has NCYCLES")
    print(f"  NCYCLES value: {plugin.container.controlParameters.NCYCLES.value}")
    print(f"  NCYCLES min: {plugin.container.controlParameters.NCYCLES.get_qualifier('min')}")


def test_prosmart_refmac_keywords_include_ncycles():
    """Test that keywordsOfTaskName includes NCYCLES."""

    # Import after setting paths
    from ccp4i2.cli.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase

    # Get keywords for prosmart_refmac
    keywords = CCP4i2RunnerBase.keywordsOfTaskName('prosmart_refmac', parent=None)

    print(f"\n=== Found {len(keywords)} keywords ===")

    # Look for NCYCLES
    ncycles_found = False
    for kw in keywords:
        if 'NCYCLES' in kw.get('path', ''):
            ncycles_found = True
            print(f"\n✓ Found NCYCLES keyword:")
            print(f"  path: {kw['path']}")
            print(f"  minimumPath: {kw['minimumPath']}")
            print(f"  qualifiers: {kw.get('qualifiers', {})}")
            break

    if not ncycles_found:
        print("\n✗ NCYCLES not found in keywords!")
        print("\nAll keyword paths:")
        for kw in keywords[:20]:  # Show first 20
            print(f"  {kw.get('path', 'NO PATH')}")
        if len(keywords) > 20:
            print(f"  ... and {len(keywords) - 20} more")

    assert ncycles_found, "NCYCLES should be in keywords list for prosmart_refmac"


if __name__ == "__main__":
    try:
        test_prosmart_refmac_plugin_has_ncycles()
        test_prosmart_refmac_keywords_include_ncycles()
        print("\n✅ All plugin keyword tests passed!")
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
