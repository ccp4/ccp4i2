"""
Test that file info persists in input_params.xml after setting cell parameters.

Reproduces the GUI workflow where:
1. Upload unmerged MTZ file to aimless_pipe UNMERGEDFILES
2. Frontend digests the file to get cell info
3. Frontend sets cell/wavelength/crystal/dataset on the list item
4. File info (baseName, relPath, dbFileId) must still be in input_params.xml

This is a regression test for a bug where set_parameter on the list item
would cause file info to be lost from input_params.xml.
"""
import xml.etree.ElementTree as ET

import pytest

from ..base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestUnmergedFilePersistence(APITestBase):
    """Test that file data survives cell parameter updates on CImportUnmerged items."""

    task_name = "aimless_pipe"

    def test_file_survives_cell_update(self, mdm2_unmerged_mtz):
        """Upload unmerged file, set cell params, verify file info persists."""

        # 1. Create project and aimless_pipe job
        self.create_project("test_unmerged_persistence")
        self.create_job()
        from ccp4i2.db.models import Job
        job_dir = str(Job.objects.get(id=self.job_id).directory)

        # 2. Upload unmerged MTZ to UNMERGEDFILES[0].file
        upload_data = self.upload_file(
            "inputData.UNMERGEDFILES[0].file", mdm2_unmerged_mtz
        )
        # Upload response is nested JSON — just check it succeeded
        assert upload_data.get('success') or upload_data.get('baseName') or \
            'baseName' in str(upload_data), f"Upload appears to have failed: {str(upload_data)[:200]}"

        # 3. Verify file info is in input_params.xml after upload
        input_params = ET.parse(f"{job_dir}/input_params.xml")
        base_names = [e.text for e in input_params.iter('baseName') if e.text]
        assert any('unmerged' in bn.lower() or 'mdm2' in bn.lower()
                    for bn in base_names), \
            f"File baseName not found after upload. baseNames: {base_names}"

        # 4. Set cell parameters (simulating what the GUI does after digest)
        cell_update = {
            "cell": {"a": 34.15, "b": 54.81, "c": 68.0,
                     "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
            "wavelength": 1.54179,
            "crystalName": "TestCrystal",
            "dataset": "TestDataset",
        }
        self.set_param("inputData.UNMERGEDFILES[0]", cell_update)

        # 5. Verify file info is STILL in input_params.xml after cell update
        input_params = ET.parse(f"{job_dir}/input_params.xml")

        # Find baseName elements inside UNMERGEDFILES (not just any baseName)
        unmerged_list = input_params.find('.//UNMERGEDFILES')
        assert unmerged_list is not None, "UNMERGEDFILES not found in input_params.xml"

        file_elem = unmerged_list.find('.//file')
        if file_elem is None:
            # Might be structured differently - check for baseName directly
            file_elem = unmerged_list

        base_name = file_elem.findtext('.//baseName') or file_elem.findtext('baseName')
        db_file_id = file_elem.findtext('.//dbFileId') or file_elem.findtext('dbFileId')

        assert base_name and len(base_name) > 0, \
            f"File baseName LOST from input_params.xml after cell update. " \
            f"UNMERGEDFILES XML:\n{ET.tostring(unmerged_list, encoding='unicode')}"

        assert db_file_id and len(db_file_id) > 0, \
            f"dbFileId LOST from input_params.xml after cell update. " \
            f"UNMERGEDFILES XML:\n{ET.tostring(unmerged_list, encoding='unicode')}"

        # 6. Also verify cell info was saved
        cell_a = unmerged_list.findtext('.//a')
        assert cell_a and float(cell_a) == pytest.approx(34.15), \
            f"Cell 'a' not saved: {cell_a}"
