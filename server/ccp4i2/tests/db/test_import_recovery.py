"""An imported file deleted from CCP4_IMPORTED_FILES by hand.

People do this --- most often to a second copy of a large unmerged file --- and
the job naming it then failed in the first program to read it, with nothing to
say why. import_recovery puts the copy back from where it was imported from when
the original still holds the same bytes, and otherwise says which file has gone,
where it came from and what to do.

Also covers the other half of the same complaint: an upload identical to a file
already in the project now uses that file instead of keeping a second copy.
"""
import hashlib
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp

from django.test import TestCase

from ...core import CCP4ErrorHandling
from ...core.base_object.error_reporting import CErrorReport
from ...db import models
from ...lib.utils.files import import_recovery
from ...lib.utils.files.import_recovery import (
    Missing, classify, report_missing_imports, restore, restore_missing_imports,
)
from ...lib.utils.files.upload_param import _LocalPathUpload, _identical_earlier_import

PAYLOAD = b"unmerged reflections, honestly" * 1000


def _md5(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()


class _ProjectOnDisk:
    """A project directory real enough for File.path to resolve, and an
    'outside' directory standing in for where the user keeps their data."""

    def __init__(self, name):
        self.root = Path(mkdtemp())
        self.outside = Path(mkdtemp())
        (self.root / "CCP4_IMPORTED_FILES").mkdir()
        self.project = models.Project.objects.create(name=name, directory=str(self.root))
        self.job = models.Job.objects.create(
            project=self.project, number="1", task_name="aimless_pipe",
            title="job 1", status=models.Job.Status.FINISHED,
        )

    def imported(self, name="hkl.mtz", data=PAYLOAD, source_name=None,
                 source_checksum=None, on_disk=True, content=1, sub_type=1):
        if on_disk:
            (self.root / "CCP4_IMPORTED_FILES" / name).write_bytes(data)
        record = models.File.objects.create(
            job=self.job, name=name, directory=models.File.Directory.IMPORT_DIR,
            type=models.FileType.objects.get_or_create(name="application/CCP4-mtz")[0],
            annotation=f"imported {name}", content=content, sub_type=sub_type,
        )
        models.FileImport.objects.create(
            file=record, name=source_name if source_name is not None else name,
            checksum=_md5(data),
            source_checksum=source_checksum if source_checksum is not None else _md5(data),
        )
        return models.File.objects.select_related("fileimport", "job").get(pk=record.pk)

    def original(self, name="hkl.mtz", data=PAYLOAD) -> Path:
        path = self.outside / name
        path.write_bytes(data)
        return path

    def destroy(self):
        rmtree(self.root, ignore_errors=True)
        rmtree(self.outside, ignore_errors=True)


class _Id:
    def __init__(self, value):
        self.value = value

    def isSet(self):
        return bool(self.value)

    def __str__(self):
        return self.value


class _FileParam:
    def __init__(self, record, path="aimless_pipe.container.inputData.UNMERGEDFILES[0].file"):
        self.dbFileId = _Id(str(record.uuid).replace("-", ""))
        self._path = path

    def object_path(self):
        return self._path


class _Container:
    """Just enough of a container: the files it names."""

    def __init__(self, *params):
        self.params = params

    def find_all_files(self):
        return list(self.params)


class ClassifyTest(TestCase):
    def setUp(self):
        import_recovery._md5_cache.clear()
        self.home = _ProjectOnDisk("recovery")

    def tearDown(self):
        self.home.destroy()

    def test_a_file_that_is_there_is_not_reported(self):
        record = self.home.imported(source_name=str(self.home.original()))
        assert classify(record) is None

    def test_an_unchanged_original_makes_it_restorable(self):
        source = self.home.original()
        record = self.home.imported(source_name=str(source), on_disk=False)
        missing = classify(record)
        assert missing.state is Missing.RESTORABLE
        assert missing.source == source
        assert str(source) in missing.describe()

    def test_restoring_puts_back_exactly_what_was_imported(self):
        source = self.home.original()
        record = self.home.imported(source_name=str(source), on_disk=False)
        assert restore(classify(record))
        assert record.path.read_bytes() == PAYLOAD
        assert classify(record) is None
        assert not list(record.path.parent.glob("*.restoring"))

    def test_restoring_survives_the_whole_directory_being_deleted(self):
        source = self.home.original()
        record = self.home.imported(source_name=str(source), on_disk=False)
        rmtree(self.home.root / "CCP4_IMPORTED_FILES")
        assert restore(classify(record))
        assert record.path.read_bytes() == PAYLOAD

    def test_an_original_that_has_gone_too_is_reported_not_restored(self):
        source = self.home.outside / "hkl.mtz"
        record = self.home.imported(source_name=str(source), on_disk=False)
        missing = classify(record)
        assert missing.state is Missing.GONE
        assert not restore(missing)
        assert not record.path.exists()

    def test_an_original_that_has_changed_is_not_copied_back(self):
        """Putting back different bytes under the imported file's name would
        silently change what the job reads."""
        source = self.home.original(data=b"reprocessed since")
        record = self.home.imported(source_name=str(source), on_disk=False)
        missing = classify(record)
        assert missing.state is Missing.CHANGED
        assert not restore(missing)

    def test_a_derived_import_names_its_source_but_is_not_rebuilt(self):
        """An MTZ split down to some columns: the source is there, but which
        columns were taken is not recorded, so rebuilding would be a guess."""
        source = self.home.original(data=b"the whole MTZ")
        record = self.home.imported(
            source_name=str(source), data=b"just F and SIGF",
            source_checksum=_md5(b"the whole MTZ"), on_disk=False,
        )
        missing = classify(record)
        assert missing.state is Missing.DERIVED
        assert "import hkl.mtz again" in missing.describe()
        assert not restore(missing)

    def test_an_upload_has_no_original_to_look_in(self):
        """A browser upload records only the file's name; a bare name must not
        be looked for relative to wherever the server happens to be running."""
        record = self.home.imported(source_name="hkl.mtz", on_disk=False)
        missing = classify(record)
        assert missing.state is Missing.UNKNOWN
        assert missing.source is None


class ContainerTest(TestCase):
    def setUp(self):
        import_recovery._md5_cache.clear()
        self.home = _ProjectOnDisk("recovery")

    def tearDown(self):
        self.home.destroy()

    def test_a_restorable_file_is_a_warning_on_its_parameter(self):
        record = self.home.imported(source_name=str(self.home.original()), on_disk=False)
        param = _FileParam(record)
        error = CErrorReport()
        # What a mustExist file reports on its own: true, and no help at all.
        error.append(klass="CDataFile", code=102, details=f"File does not exist: {record.path}",
                     name=param.object_path(), severity=CCP4ErrorHandling.SEVERITY_ERROR)
        report_missing_imports(_Container(param), error, "aimless_pipe")
        [entry] = error.getErrors()
        assert entry["code"] == import_recovery.CODE_WILL_BE_RESTORED
        assert entry["severity"] == CCP4ErrorHandling.SEVERITY_WARNING
        assert entry["name"] == param.object_path()
        assert "copied back when the job runs" in entry["details"]

    def test_a_lost_file_blocks_the_job_and_says_what_to_do(self):
        record = self.home.imported(source_name=str(self.home.outside / "gone.mtz"),
                                    on_disk=False)
        error = CErrorReport()
        report_missing_imports(_Container(_FileParam(record)), error, "aimless_pipe")
        [entry] = error.getErrors()
        assert entry["code"] == import_recovery.CODE_CANNOT_BE_RESTORED
        assert entry["severity"] == CCP4ErrorHandling.SEVERITY_ERROR
        assert "import the file again" in entry["details"]

    def test_files_that_are_present_add_nothing(self):
        record = self.home.imported(source_name=str(self.home.original()))
        error = CErrorReport()
        report_missing_imports(_Container(_FileParam(record)), error, "aimless_pipe")
        assert error.getErrors() == []

    def test_before_a_run_what_can_be_put_back_is(self):
        restorable = self.home.imported(name="a.mtz", on_disk=False,
                                        source_name=str(self.home.original("a.mtz")))
        lost = self.home.imported(name="b.mtz", on_disk=False,
                                  source_name=str(self.home.outside / "b.mtz"))
        container = _Container(_FileParam(restorable, "t.container.inputData.A"),
                               _FileParam(lost, "t.container.inputData.B"),
                               _FileParam(restorable, "t.container.inputData.C"))
        still_missing = restore_missing_imports(container)
        assert restorable.path.read_bytes() == PAYLOAD
        assert [each.file.uuid for each in still_missing] == [lost.uuid]


class IdenticalEarlierImportTest(TestCase):
    def setUp(self):
        self.home = _ProjectOnDisk("dedupe")
        self.earlier = self.home.imported()
        self.mtz = self.earlier.type

    def tearDown(self):
        self.home.destroy()

    def test_the_same_bytes_captured_the_same_way_are_reused(self):
        found = _identical_earlier_import(self.home.job, _md5(PAYLOAD), self.mtz, 1, 1)
        assert found == self.earlier

    def test_different_content_or_subtype_is_a_different_file(self):
        """Same bytes in a slot that wants them as something else: the picker
        would not have offered the earlier file, so neither does this."""
        assert _identical_earlier_import(self.home.job, _md5(PAYLOAD), self.mtz, 2, 1) is None
        assert _identical_earlier_import(self.home.job, _md5(PAYLOAD), self.mtz, 1, 5) is None

    def test_an_earlier_copy_deleted_from_disk_is_not_reused(self):
        self.earlier.path.unlink()
        assert _identical_earlier_import(self.home.job, _md5(PAYLOAD), self.mtz, 1, 1) is None

    def test_another_project_s_copy_is_not_reused(self):
        away = _ProjectOnDisk("elsewhere")
        try:
            assert _identical_earlier_import(away.job, _md5(PAYLOAD), self.mtz, 1, 1) is None
        finally:
            away.destroy()

    def test_a_blank_checksum_matches_nothing(self):
        assert _identical_earlier_import(self.home.job, "", self.mtz, 1, 1) is None


def test_an_import_by_path_remembers_the_whole_path(tmp_path):
    """FileImport.name is filled from this: a bare name is what made a deleted
    copy unrecoverable."""
    source = tmp_path / "data" / "hkl.mtz"
    source.parent.mkdir()
    source.write_bytes(b"x")
    upload = _LocalPathUpload(source)
    assert upload.name == "hkl.mtz"
    assert Path(upload.source_path) == source.resolve()
