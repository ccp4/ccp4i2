"""
Comparison tests for ServalcatConverter vs CtruncateConverter.

This test suite demonstrates that CtruncateConverter provides the exact
same API as ServalcatConverter, allowing easy A/B comparison testing.

Both converters provide:
- ipair_to_fpair(obs_file, work_directory) -> str
- ipair_to_imean(obs_file, work_directory) -> str
- to_fmean(obs_file, work_directory) -> str

Usage:
    # Swap converters by changing a single import:
    from ccp4i2.core.conversions.servalcat_converter import ServalcatConverter as Converter
    # OR
    from ccp4i2.core.conversions.ctruncate_converter import CtruncateConverter as Converter

    # Same API for both:
    output_path = Converter.to_fmean(obs_file, work_dir)
"""

import pytest
from pathlib import Path


@pytest.fixture
def demo_data_dir():
    """Path to demo data directory."""
    return Path(__file__).parent.parent / "demo_data"


@pytest.fixture
def gamma_imean_data(demo_data_dir):
    """Gamma dataset with IMEAN data."""
    gamma_path = demo_data_dir / "gamma" / "HKLOUT_unmerged.mtz"
    if not gamma_path.exists():
        pytest.skip(f"Gamma test data not found: {gamma_path}")
    return gamma_path


@pytest.fixture
def gamma_ipair_data(demo_data_dir):
    """Gamma dataset with IPAIR data."""
    gamma_path = demo_data_dir / "gamma" / "merged_intensities_Xe.mtz"
    if not gamma_path.exists():
        pytest.skip(f"Gamma IPAIR test data not found: {gamma_path}")
    return gamma_path


class TestConverterAPICompatibility:
    """Test that both converters provide the same API."""

    def test_servalcat_imean_to_fmean(self, gamma_imean_data, tmp_path):
        """Test ServalcatConverter IMEAN → FMEAN conversion."""
        from ccp4i2.core.CCP4XtalData import CObsDataFile
        from ccp4i2.core.conversions.servalcat_converter import ServalcatConverter
        import shutil

        # Create test file
        input_mtz = tmp_path / "input.mtz"
        shutil.copy2(gamma_imean_data, input_mtz)

        # Create obs file
        obs_file = CObsDataFile()
        obs_file.setFullPath(str(input_mtz))
        obs_file.setContentFlag()

        # Convert using servalcat
        work_dir = tmp_path / "servalcat_work"
        work_dir.mkdir()
        output_path = ServalcatConverter.to_fmean(obs_file, work_directory=str(work_dir))

        # Verify output
        assert Path(output_path).exists()
        print(f"✅ ServalcatConverter output: {output_path}")

    def test_ctruncate_imean_to_fmean(self, gamma_imean_data, tmp_path):
        """Test CtruncateConverter IMEAN → FMEAN conversion."""
        from ccp4i2.core.CCP4XtalData import CObsDataFile
        from ccp4i2.core.conversions.ctruncate_converter import CtruncateConverter
        import shutil

        # Create test file
        input_mtz = tmp_path / "input.mtz"
        shutil.copy2(gamma_imean_data, input_mtz)

        # Create obs file
        obs_file = CObsDataFile()
        obs_file.setFullPath(str(input_mtz))
        obs_file.setContentFlag()

        # Convert using ctruncate
        work_dir = tmp_path / "ctruncate_work"
        work_dir.mkdir()
        output_path = CtruncateConverter.to_fmean(obs_file, work_directory=str(work_dir))

        # Verify output
        assert Path(output_path).exists()
        print(f"✅ CtruncateConverter output: {output_path}")

    def test_drop_in_replacement_pattern(self, gamma_imean_data, tmp_path):
        """
        Demonstrate drop-in replacement pattern.

        This test shows how easy it is to swap between converters
        by changing a single import alias.
        """
        from ccp4i2.core.CCP4XtalData import CObsDataFile
        import shutil

        # Prepare test data
        input_mtz = tmp_path / "input.mtz"
        shutil.copy2(gamma_imean_data, input_mtz)

        obs_file = CObsDataFile()
        obs_file.setFullPath(str(input_mtz))
        obs_file.setContentFlag()

        # Test with BOTH converters using the SAME CODE
        converters = [
            ('servalcat', 'core.conversions.servalcat_converter', 'ServalcatConverter'),
            ('ctruncate', 'core.conversions.ctruncate_converter', 'CtruncateConverter'),
        ]

        results = {}
        for name, module, class_name in converters:
            # Import converter dynamically
            import importlib
            converter_module = importlib.import_module(module)
            Converter = getattr(converter_module, class_name)

            # Use SAME API for both converters
            work_dir = tmp_path / f"{name}_work"
            work_dir.mkdir()
            output_path = Converter.to_fmean(obs_file, work_directory=str(work_dir))

            # Verify output
            assert Path(output_path).exists()
            results[name] = output_path
            print(f"✅ {name}: {output_path}")

        # Both converters should produce valid output
        assert len(results) == 2
        print("\n✅ Drop-in replacement successful - both converters work with identical API")


class TestConverterComparison:
    """Compare outputs from ServalcatConverter vs CtruncateConverter."""

    @pytest.mark.skip(reason="IPAIR data needs to be verified first")
    def test_ipair_to_fpair_comparison(self, gamma_ipair_data, tmp_path):
        """Compare IPAIR → FPAIR conversion between converters."""
        from ccp4i2.core.CCP4XtalData import CObsDataFile
        from ccp4i2.core.conversions.servalcat_converter import ServalcatConverter
        from ccp4i2.core.conversions.ctruncate_converter import CtruncateConverter
        import shutil
        import gemmi

        # Prepare test data for both converters
        results = {}

        for name, Converter in [('servalcat', ServalcatConverter), ('ctruncate', CtruncateConverter)]:
            input_mtz = tmp_path / f"{name}_input.mtz"
            shutil.copy2(gamma_ipair_data, input_mtz)

            obs_file = CObsDataFile()
            obs_file.setFullPath(str(input_mtz))
            obs_file.setContentFlag()

            work_dir = tmp_path / f"{name}_work"
            work_dir.mkdir()

            output_path = Converter.ipair_to_fpair(obs_file, work_directory=str(work_dir))

            # Read output columns
            mtz = gemmi.read_mtz_file(output_path)
            column_labels = [col.label for col in mtz.columns]

            results[name] = {
                'path': output_path,
                'columns': column_labels,
                'n_reflections': len(mtz.data)
            }

        # Compare results
        print("\n=== IPAIR → FPAIR Comparison ===")
        for name, data in results.items():
            print(f"{name}:")
            print(f"  Output: {data['path']}")
            print(f"  Columns: {data['columns']}")
            print(f"  Reflections: {data['n_reflections']}")

        # Both should have FPAIR columns
        for name in ['servalcat', 'ctruncate']:
            assert 'Fplus' in results[name]['columns'] or 'F(+)' in results[name]['columns']
            assert 'Fminus' in results[name]['columns'] or 'F(-)' in results[name]['columns']


def test_api_documentation():
    """Verify both converters expose the same methods."""
    from ccp4i2.core.conversions.servalcat_converter import ServalcatConverter
    from ccp4i2.core.conversions.ctruncate_converter import CtruncateConverter

    # Both converters should have the same public methods
    required_methods = ['ipair_to_fpair', 'ipair_to_imean', 'to_fmean']

    for method_name in required_methods:
        assert hasattr(ServalcatConverter, method_name), \
            f"ServalcatConverter missing {method_name}"
        assert hasattr(CtruncateConverter, method_name), \
            f"CtruncateConverter missing {method_name}"

    print("✅ Both converters provide the same API:")
    for method in required_methods:
        print(f"  - {method}(obs_file, work_directory)")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
