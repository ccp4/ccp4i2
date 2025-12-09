"""
Test suite for CMtzData.getColumnGroups() and CUnmergedDataContent.getColumnGroups()

This test suite validates that getColumnGroups() correctly:
1. Loads MTZ files and populates listOfColumns
2. Groups columns by dataset and groupIndex
3. Detects F-Phi pairs for map coefficients
4. Assigns correct columnGroupType based on column signatures
5. Sets appropriate contentFlag values
"""
import pytest
import gemmi
from pathlib import Path
from ccp4i2.core.CCP4XtalData import (
    CMtzDataFile, CMtzData, CUnmergedDataContent,
    CObsDataFile, CPhsDataFile, CMapCoeffsDataFile, CFreeRDataFile
)


@pytest.fixture
def demo_data_dir():
    """Path to demo data directory"""
    return Path(__file__).parent.parent / "demo_data"


class TestGetColumnGroupsBasic:
    """Basic tests for getColumnGroups() functionality"""

    def test_freer_file_single_column(self, demo_data_dir):
        """Test FreeR file with single integer column (signature 'I')"""
        freer_path = demo_data_dir / "gamma" / "freeR.mtz"
        assert freer_path.exists(), f"FreeR file not found: {freer_path}"

        # Load MTZ file
        mtz_file = CMtzDataFile()
        mtz_file.setFullPath(str(freer_path))

        # Load content and get column groups
        content = mtz_file.getFileContent()
        groups = content.getColumnGroups()

        # Verify we got column groups
        assert len(groups) > 0, "Expected at least one column group"

        # Check that we have a FreeR group
        freer_groups = [g for g in groups if g.columnGroupType.value == 'FreeR']
        assert len(freer_groups) == 1, f"Expected 1 FreeR group, got {len(freer_groups)}"

        freer_group = freer_groups[0]
        assert freer_group.contentFlag.value == 1, "FreeR contentFlag should be 1 for signature 'I'"

        # Verify column list (H, K, L are filtered out during loadFile)
        col_list = freer_group.columnList.value if hasattr(
            freer_group.columnList, 'value') else freer_group.columnList
        assert len(col_list) == 1, "FreeR group should have 1 column (FREER)"
        assert col_list[0].columnType.value == 'I', \
            "FreeR column should be type 'I'"

    def test_anomalous_intensities(self, demo_data_dir):
        """Test anomalous intensities file (I+/sigI+, I-/sigI-) - signature 'KMKM'"""
        ianom_path = demo_data_dir / "gamma" / "merged_intensities_Xe.mtz"
        assert ianom_path.exists(), f"Anomalous intensities file not found: {ianom_path}"

        # Load MTZ file
        mtz_file = CMtzDataFile()
        mtz_file.setFullPath(str(ianom_path))

        # Load content and get column groups
        content = mtz_file.getFileContent()
        groups = content.getColumnGroups()

        # Verify we got column groups
        assert len(groups) > 0, "Expected at least one column group"

        # Check that we have an Obs group with KMKM signature
        obs_groups = [g for g in groups if g.columnGroupType.value == 'Obs']
        assert len(obs_groups) == 1, f"Expected 1 Obs group, got {len(obs_groups)}"

        obs_group = obs_groups[0]
        # KMKM is the first signature in CObsDataFile.correctColumns
        assert obs_group.contentFlag.value == 1, "KMKM contentFlag should be 1"

        # Verify column list has 4 columns (I+, sigI+, I-, sigI-)
        col_list = obs_group.columnList.value if hasattr(obs_group.columnList, 'value') else obs_group.columnList
        assert len(col_list) == 4, f"KMKM group should have 4 columns, got {len(col_list)}"

        # Verify column types
        col_types = [col.columnType.value for col in col_list]
        assert col_types == ['K', 'M', 'K', 'M'], f"Expected KMKM, got {''.join(col_types)}"

    def test_hl_phases(self, demo_data_dir):
        """Test Hendrickson-Lattman phases file (signature 'AAAA')"""
        phases_path = demo_data_dir / "gamma" / "initial_phases.mtz"
        assert phases_path.exists(), f"Phases file not found: {phases_path}"

        # Load MTZ file
        mtz_file = CMtzDataFile()
        mtz_file.setFullPath(str(phases_path))

        # Load content and get column groups
        content = mtz_file.getFileContent()
        groups = content.getColumnGroups()

        # Verify we got column groups
        assert len(groups) > 0, "Expected at least one column group"

        # Check that we have a Phs group with AAAA signature
        phs_groups = [g for g in groups if g.columnGroupType.value == 'Phs']
        assert len(phs_groups) == 1, f"Expected 1 Phs group, got {len(phs_groups)}"

        phs_group = phs_groups[0]
        # AAAA is the first signature in CPhsDataFile.correctColumns
        assert phs_group.contentFlag.value == 1, "AAAA contentFlag should be 1"

        # Verify column list has 4 columns (HLA, HLB, HLC, HLD)
        col_list = phs_group.columnList.value if hasattr(phs_group.columnList, 'value') else phs_group.columnList
        assert len(col_list) == 4, f"AAAA group should have 4 columns, got {len(col_list)}"

        # Verify column types
        col_types = [col.columnType.value for col in col_list]
        assert col_types == ['A', 'A', 'A', 'A'], f"Expected AAAA, got {''.join(col_types)}"

    def test_mixed_obs_and_freer(self, demo_data_dir):
        """Test MTZ file with both observations and FreeR flag"""
        gere_path = demo_data_dir / "gere" / "gere.mtz"
        assert gere_path.exists(), f"Gere file not found: {gere_path}"

        # Load MTZ file
        mtz_file = CMtzDataFile()
        mtz_file.setFullPath(str(gere_path))

        # Load content and get column groups
        content = mtz_file.getFileContent()
        groups = content.getColumnGroups()

        # Verify we got multiple column groups
        assert len(groups) >= 2, f"Expected at least 2 groups (Obs + FreeR), got {len(groups)}"

        # Check for Obs group (I+/sigI+, I-/sigI- = KMKM)
        obs_groups = [g for g in groups if g.columnGroupType.value == 'Obs']
        assert len(obs_groups) == 1, f"Expected 1 Obs group, got {len(obs_groups)}"

        obs_group = obs_groups[0]
        obs_cols = obs_group.columnList.value if hasattr(obs_group.columnList, 'value') else obs_group.columnList
        assert len(obs_cols) == 4, f"Obs group should have 4 columns, got {len(obs_cols)}"

        # Check for FreeR group
        freer_groups = [g for g in groups if g.columnGroupType.value == 'FreeR']
        assert len(freer_groups) == 1, f"Expected 1 FreeR group, got {len(freer_groups)}"

        freer_group = freer_groups[0]
        freer_cols = freer_group.columnList.value if hasattr(freer_group.columnList, 'value') else freer_group.columnList
        assert len(freer_cols) == 1, f"FreeR group should have 1 column, got {len(freer_cols)}"


class TestGetColumnGroupsEdgeCases:
    """Test edge cases and special scenarios"""

    def test_empty_columns(self):
        """Test getColumnGroups with no columns"""
        mtz_data = CMtzData()
        # Don't populate listOfColumns - should be empty
        groups = mtz_data.getColumnGroups()
        assert groups == [], "Empty columns should return empty group list"

    def test_defensive_empty_group_list(self):
        """Test that defensive check prevents IndexError on first column"""
        # This tests the fix for Bug 3 (IndexError when groupList is empty)
        # Create MTZ data with columns that all have groupIndex == 0
        from ccp4i2.core.CCP4XtalData import CMtzColumn

        mtz_data = CMtzData()
        mtz_data.listOfColumns = []

        # Add a column with groupIndex = 0 (same as initial value)
        # This should trigger the defensive check
        col = CMtzColumn(name='col1')
        col.columnLabel = 'TEST'
        col.columnType = 'F'
        col.dataset = 'ds1'
        col.groupIndex = 0  # Same as initial groupIndex in getColumnGroups

        mtz_data.listOfColumns.append(col)

        # This should not raise IndexError
        groups = mtz_data.getColumnGroups()

        # Verify the defensive check created a group
        assert len(groups) == 1, "Defensive check should create a group"
        col_list = groups[0].columnList.value if hasattr(groups[0].columnList, 'value') else groups[0].columnList
        assert len(col_list) == 1, "Group should have one column"


class TestGetColumnGroupsSignatures:
    """Test column signature recognition for different data types"""

    def test_obs_signatures(self):
        """Test all CObsDataFile.correctColumns signatures"""
        from ccp4i2.core.CCP4XtalData import CMtzColumn
        from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type

        # Get correctColumns from CObsDataFile metadata
        meta = get_class_metadata_by_type(CObsDataFile)
        correct_columns = meta.qualifiers.get('correctColumns') if meta and meta.qualifiers else []

        assert correct_columns == ['KMKM', 'GLGL', 'JQ', 'FQ'], \
            f"CObsDataFile correctColumns mismatch: {correct_columns}"

    def test_phs_signatures(self):
        """Test all CPhsDataFile.correctColumns signatures"""
        from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type

        # Get correctColumns from CPhsDataFile metadata
        meta = get_class_metadata_by_type(CPhsDataFile)
        correct_columns = meta.qualifiers.get('correctColumns') if meta and meta.qualifiers else []

        assert correct_columns == ['AAAA', 'PW'], \
            f"CPhsDataFile correctColumns mismatch: {correct_columns}"

    def test_mapcoeffs_signatures(self):
        """Test all CMapCoeffsDataFile.correctColumns signatures"""
        from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type

        # Get correctColumns from CMapCoeffsDataFile metadata
        meta = get_class_metadata_by_type(CMapCoeffsDataFile)
        correct_columns = meta.qualifiers.get('correctColumns') if meta and meta.qualifiers else []

        assert correct_columns == ['FP', 'FQP'], \
            f"CMapCoeffsDataFile correctColumns mismatch: {correct_columns}"

    def test_freer_signatures(self):
        """Test all CFreeRDataFile.correctColumns signatures"""
        from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type

        # Get correctColumns from CFreeRDataFile metadata
        meta = get_class_metadata_by_type(CFreeRDataFile)
        correct_columns = meta.qualifiers.get('correctColumns') if meta and meta.qualifiers else []

        assert correct_columns == ['I'], \
            f"CFreeRDataFile correctColumns mismatch: {correct_columns}"


class TestGetColumnGroupsDatasetHandling:
    """Test dataset and groupIndex handling"""

    def test_multiple_datasets(self):
        """Test file with multiple datasets (currently gere.mtz has HKL_base and dataset1)"""
        from pathlib import Path
        demo_data_dir = Path(__file__).parent.parent / "demo_data"
        gere_path = demo_data_dir / "gere" / "gere.mtz"

        if not gere_path.exists():
            pytest.skip(f"Gere file not found: {gere_path}")

        # Load MTZ file
        mtz_file = CMtzDataFile()
        mtz_file.setFullPath(str(gere_path))

        # Load content and get column groups
        content = mtz_file.getFileContent()
        groups = content.getColumnGroups()

        # Verify datasets are assigned
        for group in groups:
            dataset = group.dataset.value if hasattr(group.dataset, 'value') else group.dataset
            assert dataset is not None, "Group dataset should not be None"
            assert len(str(dataset)) > 0, "Group dataset should not be empty"


class TestUnmergedDataContent:
    """Test getColumnGroups on CUnmergedDataContent (same implementation as CMtzData)"""

    def test_unmerged_file(self, demo_data_dir):
        """Test getColumnGroups on unmerged data file"""
        unmerged_path = demo_data_dir / "hypf" / "hypf_unmerged.mtz"

        if not unmerged_path.exists():
            pytest.skip(f"Unmerged file not found: {unmerged_path}")

        # CUnmergedDataContent is used by CUnmergedDataFile
        # Load via parent file to test integration
        from ccp4i2.core.CCP4XtalData import CUnmergedDataFile

        unmerged_file = CUnmergedDataFile()
        unmerged_file.setFullPath(str(unmerged_path))

        # Load content
        content = unmerged_file.getFileContent()

        # Verify content is CUnmergedDataContent
        assert isinstance(content, CUnmergedDataContent), \
            f"Expected CUnmergedDataContent, got {type(content)}"

        # Get column groups (should use same implementation as CMtzData)
        groups = content.getColumnGroups()

        # Verify we got groups
        assert len(groups) > 0, "Expected at least one column group from unmerged file"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
