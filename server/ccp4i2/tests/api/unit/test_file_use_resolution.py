"""What a file named by fileIn= / fileOut= carries into the parameter it sets.

i2run resolves such a reference to the fields that identify a File to a
CDataFile. It used to give the file's identity only, not what it holds, so a
pipeline took the file with contentFlag unset and the sub-job it handed it to
refused it: phaser_simple_phil given ``--F_SIGF fileOut=import_merged[-1].OBSOUT``
failed with "got 0, requires one of IPAIR, FPAIR, IMEAN, FMEAN". The app, picking
the same file from the project, sets content, sub-type and annotation; this
holds i2run to the same.
"""
import uuid

from ccp4i2.db import models
from ccp4i2.lib.utils.files.file_use import file_dict_for_file


def _file(test_project_path, **fields):
    test_project_path.mkdir(parents=True, exist_ok=True)
    project = models.Project.objects.create(
        name="fileuse", directory=str(test_project_path / "fileuse"))
    job = models.Job.objects.create(
        uuid=uuid.uuid4(), project=project, number="1", title="import",
        task_name="import_merged", status=models.Job.Status.FINISHED)
    mtz_type, _ = models.FileType.objects.get_or_create(
        name="application/CCP4-mtz-observed")
    return models.File.objects.create(
        uuid=uuid.uuid4(), name="OBSOUT.mtz",
        directory=models.File.Directory.JOB_DIR, type=mtz_type, job=job,
        job_param_name="OBSOUT", **fields)


def test_a_resolved_file_says_what_it_holds(test_project_path):
    the_file = _file(test_project_path, content=4, sub_type=1,
                     annotation="Mean SFs from beta_blip_P3221.mtz")
    fields = file_dict_for_file(the_file)
    assert fields["contentFlag"] == 4
    assert fields["subType"] == 1
    assert fields["annotation"] == "Mean SFs from beta_blip_P3221.mtz"
    assert fields["dbFileId"] == str(the_file.uuid).replace("-", "")


def test_what_the_record_does_not_know_is_left_out(test_project_path):
    the_file = _file(test_project_path)
    fields = file_dict_for_file(the_file)
    # Absent, not None: a None would unset what the file itself declares.
    assert "contentFlag" not in fields and "subType" not in fields


def test_a_bare_db_file_id_resolves_to_the_same_fields(test_project_path):
    # i2run's dbFileId= gave the id alone, so the file had no path when the
    # job was validated and its content was never checked: a PDB-format
    # model passed for a task that needs mmCIF, and the job failed later.
    from ccp4i2.lib.utils.files.file_use import resolve_db_file_id

    the_file = _file(test_project_path, content=4, sub_type=1)
    assert resolve_db_file_id(str(the_file.uuid)) == file_dict_for_file(the_file)


def test_an_unknown_db_file_id_is_an_error(test_project_path):
    import pytest

    from ccp4i2.lib.utils.files.file_use import FileUseError, resolve_db_file_id

    _file(test_project_path)
    with pytest.raises(FileUseError):
        resolve_db_file_id(str(uuid.uuid4()))


def test_resolve_fileuse_takes_a_file_id(test_project_path):
    # A file that went into a list element has no [job].PARAM reference that
    # resolves; its id names it, so an agent can give it back.
    from ccp4i2.lib.utils.files.resolve_fileuse import resolve_fileuse
    the_file = _file(test_project_path, content=4)
    project = the_file.job.project
    for text in (str(the_file.uuid), str(the_file.uuid).replace("-", "")):
        result = resolve_fileuse(project, text)
        assert result.success, result.error
        assert result.data["dbFileId"] == str(the_file.uuid).replace("-", "")
        assert result.data["baseName"] == "OBSOUT.mtz"
    other = models.Project.objects.create(name="other", directory=str(test_project_path / "other"))
    assert not resolve_fileuse(other, str(the_file.uuid)).success   # not across projects


def test_a_consumed_name_resolves_to_the_jobs_out_file_when_it_wrote_one(test_project_path):
    # Haiku followed the judgement's "[6].F_SIGF" into the original P 21 2 21
    # file after Phaser had changed the group to P 21 21 21 and written
    # F_SIGF_OUT (Eleanor's case, 2026-10-10). Downstream of such a job the
    # consumed name means the job's own re-expression of it.
    from ccp4i2.lib.utils.files.resolve_fileuse import resolve_fileuse
    original = _file(test_project_path, content=4, sub_type=1, annotation="/BAD1330/NATIVE")
    project = original.job.project
    phaser = models.Job.objects.create(
        uuid=uuid.uuid4(), project=project, number="2", title="Phaser basic",
        task_name="phaser_simple_phil", status=models.Job.Status.FINISHED)
    models.FileUse.objects.create(file=original, job=phaser, role=models.FileUse.Role.IN,
                                  job_param_name="F_SIGF")
    relabelled = models.File.objects.create(
        uuid=uuid.uuid4(), name="F_SIGF_OUT.mtz", directory=models.File.Directory.JOB_DIR,
        type=original.type, job=phaser, job_param_name="F_SIGF_OUT", content=4, sub_type=1,
        annotation="/BAD1330/NATIVE, not reindexed (P 21 21 21, 1.26 A)")
    assert resolve_fileuse(project, "[2].F_SIGF").data["baseName"] == "F_SIGF_OUT.mtz"
    assert resolve_fileuse(project, "phaser_simple_phil[-1].F_SIGF").data["dbFileId"] == relabelled.uuid.hex
    # explicit names still mean what they say, and a job without an _OUT is untouched
    assert resolve_fileuse(project, "[2].F_SIGF_OUT").data["baseName"] == "F_SIGF_OUT.mtz"
    assert resolve_fileuse(project, "[1].OBSOUT").data["baseName"] == "OBSOUT.mtz"
