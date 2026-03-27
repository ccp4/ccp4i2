"""
Test that qualifiers are properly inherited through the CData class hierarchy.

This test verifies that fileContentClassName set in CMtzDataFile.__init__()
is accessible in subclasses like CObsDataFile, CMiniMtzDataFile, etc.
"""


def test_qualifier_inheritance():
    """Test that fileContentClassName is inherited from CMtzDataFile."""

    from ccp4i2.core.CCP4XtalData import (
        CMtzDataFile,
        CMiniMtzDataFile,
        CObsDataFile,
        CPhsDataFile,
        CFreeRDataFile
    )

    # Test CMtzDataFile (base class that sets the qualifier)
    mtz = CMtzDataFile()
    content_class = mtz.get_qualifier('fileContentClassName')
    print(f"CMtzDataFile fileContentClassName: {content_class}")
    assert content_class == 'CMtzData', f"Expected 'CMtzData', got {content_class}"

    # Test CMiniMtzDataFile (direct subclass)
    mini_mtz = CMiniMtzDataFile()
    content_class = mini_mtz.get_qualifier('fileContentClassName')
    print(f"CMiniMtzDataFile fileContentClassName: {content_class}")
    assert content_class == 'CMtzData', f"Expected 'CMtzData', got {content_class}"

    # Test CObsDataFile (subclass of CMiniMtzDataFile)
    obs = CObsDataFile()
    content_class = obs.get_qualifier('fileContentClassName')
    print(f"CObsDataFile fileContentClassName: {content_class}")
    assert content_class == 'CMtzData', f"Expected 'CMtzData', got {content_class}"

    # Test CPhsDataFile (subclass of CMiniMtzDataFile)
    phs = CPhsDataFile()
    content_class = phs.get_qualifier('fileContentClassName')
    print(f"CPhsDataFile fileContentClassName: {content_class}")
    assert content_class == 'CMtzData', f"Expected 'CMtzData', got {content_class}"

    # Test CFreeRDataFile (subclass of CMiniMtzDataFile)
    free_r = CFreeRDataFile()
    content_class = free_r.get_qualifier('fileContentClassName')
    print(f"CFreeRDataFile fileContentClassName: {content_class}")
    assert content_class == 'CMtzData', f"Expected 'CMtzData', got {content_class}"

    print("\n✅ All subclasses correctly inherit fileContentClassName='CMtzData'")


if __name__ == "__main__":
    try:
        test_qualifier_inheritance()
        print("\n🎉 Qualifier inheritance test passed!")
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
