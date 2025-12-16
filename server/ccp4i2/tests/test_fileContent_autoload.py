"""Test that fileContent auto-loads when accessing the property."""

from pathlib import Path

from ccp4i2.core import CCP4Utils


def test_fileContent_autoload_with_valid_path():
    """Test that fileContent auto-loads when file has a valid path."""
    from ccp4i2.core.CCP4XtalData import CObsDataFile

    # Use the demo data file
    demo_file = Path(CCP4Utils.getCCP4I2Dir()) / "demo_data" / "gamma" / "merged_intensities_Xe.mtz"
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
