"""Test that fileContent auto-loads when accessing the property."""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "server"))

# Set CCP4I2_ROOT for plugin discovery
os.environ["CCP4I2_ROOT"] = str(project_root)


def test_fileContent_autoload_with_valid_path():
    """Test that fileContent auto-loads when file has a valid path."""
    from core.CCP4XtalData import CObsDataFile
    from pathlib import Path

    # Use the demo data file
    demo_file = Path(project_root) / "demo_data" / "gamma" / "merged_intensities_Xe.mtz"
    assert demo_file.exists(), f"Demo file not found: {demo_file}"

    # Create file object and set path
    obs_file = CObsDataFile()
    obs_file.baseName = str(demo_file)

    # Access fileContent property (should auto-create and auto-load)
    content = obs_file.fileContent

    assert content is not None, "fileContent should be auto-created"
    assert content.resolutionRange is not None, "resolutionRange should be initialized"
    assert content.resolutionRange.high is not None, "resolutionRange.high should be loaded"

    high_res = float(content.resolutionRange.high)
    print(f"High resolution: {high_res}")
    assert high_res > 0, f"High resolution should be > 0, got {high_res}"

    print(f"✅ Auto-load worked! Resolution range: {content.resolutionRange.low} - {content.resolutionRange.high} Å")


if __name__ == "__main__":
    try:
        test_fileContent_autoload_with_valid_path()
        print("\n🎉 fileContent auto-load test passed!")
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
