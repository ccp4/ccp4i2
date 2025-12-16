"""
Test that prosmart_refmac accepts --NCYCLES argument via i2run.

This is the end-to-end test to verify the hash collision fix works
for actual command-line argument parsing.
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

from ccp4i2.cli.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase


def test_prosmart_refmac_has_ncycles_argument():
    """Test that prosmart_refmac plugin accepts --NCYCLES argument."""

    # Get keywords for prosmart_refmac
    keywords = CCP4i2RunnerBase.keywordsOfTaskName('prosmart_refmac', parent=None)

    print(f"Found {len(keywords)} keywords for prosmart_refmac")

    # Look for NCYCLES
    ncycles_found = False
    ncycles_keyword = None
    for kw in keywords:
        if 'NCYCLES' in kw.get('path', ''):
            ncycles_found = True
            ncycles_keyword = kw
            break

    if ncycles_found:
        print(f"\n✅ NCYCLES found in keywords!")
        print(f"  path: {ncycles_keyword['path']}")
        print(f"  minimumPath: {ncycles_keyword.get('minimumPath', 'N/A')}")
        print(f"  className: {ncycles_keyword.get('className', 'N/A')}")
        print(f"  qualifiers: {ncycles_keyword.get('qualifiers', {})}")
    else:
        print(f"\n❌ NCYCLES not found in keywords!")
        print(f"\nAll keyword paths:")
        for kw in keywords[:30]:
            print(f"  {kw.get('path', 'NO PATH')}")
        if len(keywords) > 30:
            print(f"  ... and {len(keywords) - 30} more")

    assert ncycles_found, "NCYCLES should be in keywords list for prosmart_refmac"


def test_all_missing_parameters_present():
    """Test that all 14 previously missing parameters are now present."""

    # These were the 14 parameters that were missing due to hash collision
    expected_params = [
        'NCYCLES',
        'NTLSCYCLES_AUTO',
        'PROSMART_PROTEIN_WEIGHT',
        'PROSMART_PROTEIN_ALPHA',
        'PROSMART_PROTEIN_DMAX',
        'PROSMART_NUCLEICACID_WEIGHT',
        'PROSMART_NUCLEICACID_ALPHA',
        'PROSMART_NUCLEICACID_DMAX',
        'BSHARP',
        'RESOLUTION',
        'RES_MIN',
        'RES_MAX',
        'WAVELENGTH',
        'SOLVENT_SHRINK'
    ]

    # Get keywords for prosmart_refmac
    keywords = CCP4i2RunnerBase.keywordsOfTaskName('prosmart_refmac', parent=None)

    # Build set of parameter names from keyword paths
    param_names = set()
    for kw in keywords:
        path = kw.get('path', '')
        # Extract parameter name from path (e.g., "prosmart_refmac.controlParameters.NCYCLES" -> "NCYCLES")
        if '.' in path:
            param_name = path.split('.')[-1]
            param_names.add(param_name)

    print(f"\n=== Checking for previously missing parameters ===")
    missing = []
    found = []

    for param in expected_params:
        if param in param_names:
            found.append(param)
            print(f"  ✅ {param}")
        else:
            missing.append(param)
            print(f"  ❌ {param}")

    print(f"\nSummary: {len(found)}/14 found, {len(missing)}/14 still missing")

    if missing:
        print(f"\nStill missing: {', '.join(missing)}")

    assert len(missing) == 0, f"All 14 parameters should be present, but {len(missing)} are still missing: {missing}"


if __name__ == "__main__":
    try:
        test_prosmart_refmac_has_ncycles_argument()
        test_all_missing_parameters_present()
        print("\n🎉 All i2run NCYCLES tests passed!")
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
