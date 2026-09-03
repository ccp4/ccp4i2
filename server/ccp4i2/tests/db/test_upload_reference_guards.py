"""Nothing may delete a file another job still depends on.

The upload path replaces the file a parameter holds, and replacing means
unlinking the old one. `_is_referenced_elsewhere` guards that, but it only sees
the container being edited --- so it cannot see a *second job* pointing at the
same file. That happens routinely: the frontend's "Browse the project
hierarchy" dialog sets dbFileId to a file owned by another job, in this project
or another one. Replacing the owning job's upload then unlinked bytes that job
was still using.

`_is_referenced_outside_this_job` answers from two sources, because neither is
complete on its own: FileUse rows for jobs that have run, and the params XML of
jobs that have not (which have no FileUse row yet).

Also covers the project scoping of duplicate detection. File.path resolves
through file.job.project.directory, so a match from another project would point
this job at bytes under a project it does not own.
"""
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp

from django.test import TestCase

from ...db import models
from ...lib.utils.files.upload_param import (
    _discard_staged_upload, _existing_imports_of_same_source,
    _is_referenced_outside_this_job,
)


class _ProjectOnDisk:
    """A project directory real enough for File.path to resolve."""

    def __init__(self, name):
        self.root = Path(mkdtemp())
        (self.root / "CCP4_IMPORTED_FILES").mkdir()
        self.project = models.Project.objects.create(
            name=name, directory=str(self.root)
        )

    def job(self, number, status=models.Job.Status.FINISHED):
        job = models.Job.objects.create(
            project=self.project, number=number, task_name="prosmart_refmac",
            title=f"job {number}", status=status,
        )
        job.directory.mkdir(parents=True, exist_ok=True)
        return job

    def imported_file(self, job, name, source_checksum="", on_disk=True,
                      checksum="0" * 32):
        if on_disk:
            (self.root / "CCP4_IMPORTED_FILES" / name).write_bytes(b"content")
        file_record = models.File.objects.create(
            job=job, name=name, directory=models.File.Directory.IMPORT_DIR,
            type=models.FileType.objects.get_or_create(name="application/CCP4-mtz")[0],
            annotation=f"imported {name}",
        )
        models.FileImport.objects.create(
            file=file_record, name=name, checksum=checksum,
            source_checksum=source_checksum,
        )
        return file_record

    def destroy(self):
        rmtree(self.root, ignore_errors=True)


class ReferenceGuardTest(TestCase):
    def setUp(self):
        self.home = _ProjectOnDisk("guards")
        self.owner = self.home.job("1")
        self.file = self.home.imported_file(self.owner, "gamma.mtz")

    def tearDown(self):
        self.home.destroy()

    def test_a_file_nobody_else_uses_may_be_deleted(self):
        assert not _is_referenced_outside_this_job(self.file, self.owner)

    def test_the_owning_job_s_own_use_does_not_count(self):
        """Otherwise a job could never replace a file it had already run with."""
        models.FileUse.objects.create(
            file=self.file, job=self.owner, role=models.FileUse.Role.IN,
            job_param_name="F_SIGF",
        )
        assert not _is_referenced_outside_this_job(self.file, self.owner)

    def test_a_file_another_job_has_run_with_is_kept(self):
        consumer = self.home.job("2")
        models.FileUse.objects.create(
            file=self.file, job=consumer, role=models.FileUse.Role.IN,
            job_param_name="F_SIGF",
        )
        assert _is_referenced_outside_this_job(self.file, self.owner)

    def test_a_job_in_another_project_counts_too(self):
        """Selected through "Browse the project hierarchy", which is not
        restricted to the current project."""
        away = _ProjectOnDisk("elsewhere")
        try:
            models.FileUse.objects.create(
                file=self.file, job=away.job("1"), role=models.FileUse.Role.IN,
                job_param_name="F_SIGF",
            )
            assert _is_referenced_outside_this_job(self.file, self.owner)
        finally:
            away.destroy()

    def test_a_pending_job_holding_the_reference_in_its_params_is_seen(self):
        """The case FileUse cannot cover: the job has not run, so no FileUse
        row exists, but its params already name the file."""
        pending = self.home.job("2", status=models.Job.Status.PENDING)
        (pending.directory / "input_params.xml").write_text(
            f"<params><F_SIGF><dbFileId>{self.file.uuid}</dbFileId></F_SIGF></params>"
        )
        assert _is_referenced_outside_this_job(self.file, self.owner)

    def test_the_undashed_uuid_form_is_seen_too(self):
        """What actually gets written: dbFileId is stored with the dashes
        stripped, so matching only the dashed form would miss every real case."""
        pending = self.home.job("2", status=models.Job.Status.PENDING)
        undashed = str(self.file.uuid).replace("-", "")
        (pending.directory / "input_params.xml").write_text(
            f"<params><F_SIGF><dbFileId>{undashed}</dbFileId></F_SIGF></params>"
        )
        assert _is_referenced_outside_this_job(self.file, self.owner)

    def test_a_pending_job_naming_some_other_file_does_not_block(self):
        pending = self.home.job("2", status=models.Job.Status.PENDING)
        (pending.directory / "input_params.xml").write_text(
            "<params><F_SIGF><dbFileId>"
            "0123456789abcdef0123456789abcdef"
            "</dbFileId></F_SIGF></params>"
        )
        assert not _is_referenced_outside_this_job(self.file, self.owner)

    def test_a_finished_job_s_params_are_not_scanned(self):
        """Jobs that have run are answered for by FileUse. Scanning their
        params as well would be redundant work on every upload."""
        ran = self.home.job("2", status=models.Job.Status.FINISHED)
        (ran.directory / "input_params.xml").write_text(
            f"<params><dbFileId>{self.file.uuid}</dbFileId></params>"
        )
        assert not _is_referenced_outside_this_job(self.file, self.owner)

    def test_an_unanswerable_question_keeps_the_file(self):
        """False means "go ahead and unlink", which is the wrong way for an
        unanswerable question about a destructive act to fail.

        An unsaved record cannot be used in a related filter --- Django raises
        rather than returning nothing --- which stands in here for any reason
        the reference check might not complete.
        """
        unsaved = models.File(name="x", directory=models.File.Directory.IMPORT_DIR)
        assert _is_referenced_outside_this_job(unsaved, self.owner)


class DuplicateAdvisoryTest(TestCase):
    def setUp(self):
        self.home = _ProjectOnDisk("advisory")
        self.job = self.home.job("1")

    def tearDown(self):
        self.home.destroy()

    def test_an_earlier_import_of_the_same_bytes_is_reported(self):
        earlier = self.home.imported_file(self.job, "gamma.mtz", source_checksum="abc")
        found = _existing_imports_of_same_source(self.job, "abc")
        assert [each["uuid"] for each in found] == [str(earlier.uuid)]
        assert found[0]["baseName"] == "gamma.mtz"
        assert found[0]["jobNumber"] == "1"

    def test_different_bytes_are_not_reported(self):
        self.home.imported_file(self.job, "gamma.mtz", source_checksum="abc")
        assert _existing_imports_of_same_source(self.job, "def") == []

    def test_a_match_in_another_project_is_not_offered(self):
        """File.path resolves through file.job.project.directory, so reusing
        this would point the job at bytes under a project it does not own."""
        away = _ProjectOnDisk("away")
        try:
            away.imported_file(away.job("1"), "gamma.mtz", source_checksum="abc")
            assert _existing_imports_of_same_source(self.job, "abc") == []
        finally:
            away.destroy()

    def test_a_match_whose_file_has_gone_from_disk_is_not_offered(self):
        self.home.imported_file(
            self.job, "gamma.mtz", source_checksum="abc", on_disk=False
        )
        assert _existing_imports_of_same_source(self.job, "abc") == []

    def test_the_file_just_written_can_be_excluded(self):
        """The upload records its own FileImport before the advisory is
        returned, and must not report the upload as a duplicate of itself."""
        just_written = self.home.imported_file(
            self.job, "gamma.mtz", source_checksum="abc"
        )
        assert _existing_imports_of_same_source(
            self.job, "abc", exclude_file=just_written
        ) == []


class SourceChecksumRoundTripTest(TestCase):
    """A project that has been exported and re-imported must not silently lose
    the ability to spot re-uploads of files it already holds."""

    def setUp(self):
        self.home = _ProjectOnDisk("roundtrip")
        self.job = self.home.job("1")

    def tearDown(self):
        self.home.destroy()

    def _exported_element(self):
        import xml.etree.ElementTree as ET

        from ...db.export_project import _export_file_import_table

        body = ET.Element("body")
        _export_file_import_table(body, self.home.project)
        return body.find("importfileTable/importfile")

    def test_the_source_checksum_survives_export_and_import(self):
        from ...db.import_i2xml import import_file_import

        self.home.imported_file(self.job, "gamma.mtz", source_checksum="abc123")
        element = self._exported_element()
        assert element.get("sourcechecksum") == "abc123"

        models.FileImport.objects.all().delete()
        import_file_import(element)
        assert models.FileImport.objects.get().source_checksum == "abc123"

    def test_an_archive_written_before_the_field_existed_still_imports(self):
        """Legacy CCP4i2 exports, and our own from before this change, carry no
        sourcechecksum attribute. Blank means "takes no part in duplicate
        detection", which is the right answer for them."""
        from ...db.import_i2xml import import_file_import

        self.home.imported_file(self.job, "gamma.mtz", source_checksum="abc123")
        element = self._exported_element()
        del element.attrib["sourcechecksum"]

        models.FileImport.objects.all().delete()
        import_file_import(element)
        assert models.FileImport.objects.get().source_checksum == ""


class InterchangeabilityTest(TestCase):
    """Same source bytes does not mean interchangeable.

    Reconstructed from project 1h1s_new. The PDB reflection file r1h1ssf.ent
    was fetched once for F/SIGF and again for the free-R set. Both fetches
    staged byte-identical sources, but the two imported files hold different
    columns:

        CCP4_DOWNLOADED_FILES  r1h1ssf.ent    037b9fb1...  }  identical
                               r1h1ssf_1.ent  037b9fb1...  }
        CCP4_IMPORTED_FILES    r1h1ssf.mtz    0c5e22be...  FP, SIGFP
                               r1h1ssf_1.mtz  a9581d7f...  FreeR_flag

    Under the current design that second fetch is the only way to get the
    free-R set --- each import derives one representation and the source is
    not kept as something a second derivation could be taken from. Reporting
    it as a redundant upload would be wrong, and the fastest way to teach
    people to ignore the message.
    """

    def setUp(self):
        self.home = _ProjectOnDisk("interchangeability")
        self.job = self.home.job("1")

    def tearDown(self):
        self.home.destroy()

    def _typed(self, name, type_name, content, source_checksum):
        file_record = models.File.objects.create(
            job=self.job, name=name, directory=models.File.Directory.IMPORT_DIR,
            type=models.FileType.objects.get_or_create(name=type_name)[0],
            content=content, annotation=name,
        )
        (self.home.root / "CCP4_IMPORTED_FILES" / name).write_bytes(b"x")
        models.FileImport.objects.create(
            file=file_record, name="r1h1ssf.ent", checksum="0" * 32,
            source_checksum=source_checksum,
        )
        return file_record

    def test_both_representations_are_found_as_same_source(self):
        """The detection itself is right --- they *are* the same bytes."""
        self._typed("r1h1ssf.mtz", "application/CCP4-mtz", 4, "037b9fb1")
        found = _existing_imports_of_same_source(self.job, "037b9fb1")
        assert [each["baseName"] for each in found] == ["r1h1ssf.mtz"]

    def test_a_different_representation_is_reported_but_not_as_reusable(self):
        """What the free-R fetch must produce: seen, not scolded."""
        self._typed("r1h1ssf.mtz", "application/CCP4-mtz", 4, "037b9fb1")
        found = _existing_imports_of_same_source(self.job, "037b9fb1")
        # The classification the upload path applies against the new file:
        # same type, different content flag.
        for each in found:
            each["interchangeable"] = (
                each["type"] == "application/CCP4-mtz" and each["contentFlag"] == 0
            )
        assert found and not any(each["interchangeable"] for each in found)

    def test_the_same_representation_is_reusable(self):
        """The case the advisory exists for: this really is a re-upload."""
        self._typed("gamma.mtz", "application/CCP4-mtz", 4, "abc123")
        found = _existing_imports_of_same_source(self.job, "abc123")
        for each in found:
            each["interchangeable"] = (
                each["type"] == "application/CCP4-mtz" and each["contentFlag"] == 4
            )
        assert found and all(each["interchangeable"] for each in found)


class StagedUploadCleanupTest(TestCase):
    """A failed import must not leave its bytes in the project.

    Reflection uploads stage into CCP4_DOWNLOADED_FILES before anything is
    derived from them, and plenty of things can fail after that: an unmerged
    file offered to a parameter needing merged data, an MTZ with no column
    group matching the requested type. Until a File row was written the staged
    bytes had no owner, so a failure left them with nothing referencing them
    and nothing that would ever collect them --- and nothing cleans that
    directory, ever.

    Measured: two failed uploads of one 5.9 MB unmerged file left 11.8 MB.
    """

    def setUp(self):
        self.home = _ProjectOnDisk("staging")
        self.job = self.home.job("1")
        self.staged_dir = self.home.root / "CCP4_DOWNLOADED_FILES"
        self.staged_dir.mkdir(exist_ok=True)

    def tearDown(self):
        self.home.destroy()

    def _stage(self, name="gamma_native.mtz"):
        path = self.staged_dir / name
        path.write_bytes(b"MTZ staged bytes")
        return path

    def test_a_staged_file_nothing_refers_to_is_removed(self):
        staged = self._stage()
        _discard_staged_upload(self.job, staged)
        assert not staged.exists()

    def test_a_staged_file_a_File_row_names_is_kept(self):
        """A failure *after* the row was written must not delete the file the
        database now points at."""
        staged = self._stage("gamma.mtz")
        self.home.imported_file(self.job, "gamma.mtz")
        _discard_staged_upload(self.job, staged)
        assert staged.exists()

    def test_a_row_in_another_project_does_not_protect_it(self):
        """Scoped to this project: an unrelated project holding a file of the
        same name says nothing about these bytes."""
        away = _ProjectOnDisk("staging_away")
        try:
            staged = self._stage("gamma.mtz")
            away.imported_file(away.job("1"), "gamma.mtz")
            _discard_staged_upload(self.job, staged)
            assert not staged.exists()
        finally:
            away.destroy()

    def test_a_file_already_gone_is_not_an_error(self):
        """The no-split path moves the staged file rather than copying it, so
        by the time a later step fails it may not be there at all."""
        _discard_staged_upload(self.job, self.staged_dir / "never-existed.mtz")

    def test_cleanup_never_raises(self):
        """A failed cleanup must never replace the error that caused it."""
        _discard_staged_upload(self.job, None)
        _discard_staged_upload(self.job, self.staged_dir)  # a directory
