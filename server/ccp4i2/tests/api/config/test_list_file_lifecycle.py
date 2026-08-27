"""
Configuring a job through the API, without running it.

A third API tier. `api/unit` checks endpoint behaviour and `api/e2e` runs whole
pipelines; neither exercises what a user actually does before pressing Run ---
adding rows to a list, uploading into them, removing one, adding another. That
gap is why job 66 of GammaBySAD lost a file: the sequence needs add, add,
remove, re-add through the REST API, and nothing did that.

No task execution, so no CCP4 binaries and no minutes: these are seconds.

The invariant under test is simple and was violated:

    replacing the file in one slot must never delete a file another slot
    still refers to
"""
import pytest

from ..base import APITestBase


class TestListFileLifecycle(APITestBase):
    """Rows come and go; the files behind other rows must not."""

    task_name = "aimless_pipe"

    def _files_on_disk(self):
        return sorted(p.name for p in
                      (self.get_job_directory().parent.parent / "CCP4_IMPORTED_FILES").glob("*.mtz"))

    def test_replacing_one_row_leaves_another_rows_file_alone(
            self, gamma_unmerged_mtz, mdm2_unmerged_mtz):
        """The shape of the job 66 failure.

        Slot [1]'s record was found by the positional name 'UNMERGEDFILES[1]',
        which after a removal belonged to a different row --- so uploading into
        [1] unlinked the file that [0] was still pointing at.
        """
        self.create_project("config_replace_row")
        self.create_job()

        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", gamma_unmerged_mtz)
        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[1].file", mdm2_unmerged_mtz)

        before = self._files_on_disk()
        assert len(before) == 2, f"expected two imported files, got {before}"

        # Replace what is in the second slot.
        self.upload_file("inputData.UNMERGEDFILES[1].file", mdm2_unmerged_mtz)

        after = self._files_on_disk()
        assert before[0] in after or any(n.startswith('gamma') for n in after), \
            f"the first row's file was deleted by an upload into the second: {after}"

    def test_the_same_file_in_two_rows_survives_one_being_replaced(
            self, gamma_unmerged_mtz, mdm2_unmerged_mtz):
        """A user may legitimately add the same data twice."""
        self.create_project("config_same_file_twice")
        self.create_job()

        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", gamma_unmerged_mtz)
        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[1].file", gamma_unmerged_mtz)

        self.upload_file("inputData.UNMERGEDFILES[1].file", mdm2_unmerged_mtz)

        remaining = self._files_on_disk()
        assert any(n.startswith('gamma') for n in remaining), \
            f"the first row lost its file when the second was replaced: {remaining}"

    def _list_value(self, name="UNMERGEDFILES"):
        """The list as the client sees it, so a row can be spliced out.

        The container endpoint nests every node as {_class, _value}; the
        client reads `item._value` for a list, which is the array of rows.
        """
        response = self.client.get(f'{self.API_PREFIX}/jobs/{self.job_id}/container/')
        assert response.status_code == 200, response.content
        node = response.json()["data"]["result"]["_value"]["inputData"]["_value"][name]
        return node["_value"]

    def test_the_reported_sequence_add_add_remove_readd(
            self, gamma_unmerged_mtz, mdm2_unmerged_mtz):
        """Job 66 of GammaBySAD, reproduced.

        Add gamma, add MDM2, remove a row, add MDM2 back. The gamma row then
        read "File does not exist", its dbFileId matched no record, and its
        file was gone from CCP4_IMPORTED_FILES.
        """
        from ccp4i2.db import models

        self.create_project("config_add_add_remove_readd")
        self.create_job()

        # Three rows, and the middle one goes. That is what makes the records
        # and the rows disagree: removing the last renumbers nothing.
        for index, path in enumerate((gamma_unmerged_mtz, mdm2_unmerged_mtz,
                                      gamma_unmerged_mtz)):
            self.add_list_item("inputData.UNMERGEDFILES")
            self.upload_file(f"inputData.UNMERGEDFILES[{index}].file", path)

        records = {f.job_param_name: f for f in
                   models.File.objects.filter(job__id=self.job_id)}
        assert len(records) == 3, f"expected three imports, got {sorted(records)}"
        last = records["UNMERGEDFILES[2].file"]
        survivor_path = last.path
        survivor_uuid = last.uuid

        current = self._list_value()
        assert len(current) == 3, f"expected three rows, got {len(current)}"
        self.set_param("inputData.UNMERGEDFILES", [current[0], current[2]])
        assert len(self._list_value()) == 2, \
            "the row was not removed, so this never reaches the case under test"

        # The third row now sits at index 1, while its record still says [2].
        # Adding another row appends at index 2 and uploads into it.
        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[2].file", mdm2_unmerged_mtz)

        gamma_path = survivor_path
        gamma_records = [last]
        assert models.File.objects.filter(uuid=survivor_uuid).exists(), \
            "a still-referenced row's record was deleted by an upload into a " \
            "slot that merely shares its positional name"

        assert gamma_path.exists(), (
            "the gamma row's file was deleted by re-adding MDM2 --- a file "
            "record was matched by its positional name rather than its identity")
        assert models.File.objects.filter(uuid=gamma_records[0].uuid).exists(), \
            "the gamma row now points at a record that no longer exists"


    def test_every_row_still_points_at_a_file_that_exists(
            self, gamma_unmerged_mtz, mdm2_unmerged_mtz):
        """The symptom the user saw: a row referring to a deleted record.

        'File does not exist' with a dbFileId matching nothing is the visible
        end of the same fault; this asserts the state rather than the message.
        """
        self.create_project("config_rows_resolve")
        self.create_job()

        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", gamma_unmerged_mtz)
        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[1].file", mdm2_unmerged_mtz)
        self.upload_file("inputData.UNMERGEDFILES[1].file", mdm2_unmerged_mtz)

        imported = self.get_job_directory().parent.parent / "CCP4_IMPORTED_FILES"
        from ccp4i2.db import models

        for record in models.File.objects.filter(job__id=self.job_id):
            assert (imported / record.name).exists() or record.path.exists(), \
                f"{record.name} is in the database but not on disk"
