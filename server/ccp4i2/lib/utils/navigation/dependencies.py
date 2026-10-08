"""Which jobs go when a job is deleted, and deleting them.

Deleting a job deletes everything that depends on it: its sub-jobs, and every
job that used one of its files, recursively.

Files a job *imported* (a FileImport row each, the bytes in the project's
CCP4_IMPORTED_FILES) are the one choice the user gets, as in Qt-i2: delete
them with the job, or keep them. Keeping is the default. A kept file stays in
the project and stays usable, so the jobs that used it are not dependents and
survive. File.job cascades, and an imported file's path is resolved through
its job's project, so its row cannot outlive every job: like Qt-i2
(CCP4DbApi.setJobToImport), the importing job is kept as a FILE_HOLDER --- its
outputs, sub-jobs and working files gone, its imported files still attached.
"""

from typing import Callable, Dict, Iterable, List, Optional, Set
import logging
import shutil
from pathlib import Path

from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")
logger.setLevel(logging.WARNING)

ACTIVE_STATUSES = {
    models.Job.Status.QUEUED,
    models.Job.Status.RUNNING,
    models.Job.Status.RUNNING_REMOTELY,
}

#: What a file holder keeps of its job directory: the parameter files, so the
#: job still opens and its input still names the files it holds.
FILE_HOLDER_KEEPS = ("params.xml", "input_params.xml")


# Note that python handles tuple sort exactly as one would wish
# https://stackoverflow.com/questions/2574080/how-to-sort-in-python-a-list-of-strings-containing-numbers
def version_sort_key(job: models.Job) -> tuple:
    return tuple(map(int, job.number.split(".")))


def imported_files_of(job: models.Job):
    """The files ``job`` brought into the project (a FileImport row each).

    They live in CCP4_IMPORTED_FILES, not in the job's directory, and other
    jobs may use them: an upload whose bytes match an earlier import reuses
    that File row rather than copying it again.
    """
    return models.File.objects.filter(job=job, fileimport__isnull=False)


def _follows_imported_files(job: models.Job, follow_imported_files: bool) -> bool:
    # A file holder exists only to hold its imported files, so deleting one
    # always means deleting them; keeping them would leave it where it was.
    return follow_imported_files or job.status == models.Job.Status.FILE_HOLDER


def find_dependent_jobs(
    the_job: models.Job,
    growing_list: List[models.Job] = None,
    leaf_action=None,
    follow_imported_files: bool = True,
) -> List[models.Job]:
    """Jobs that must go if ``the_job`` goes: its sub-jobs, and every job that
    used one of its files, recursively.

    ``follow_imported_files`` decides whether a job that merely used a file
    ``the_job`` imported counts. It does when the imported files are deleted
    with the job; it does not when they are kept, since the file stays usable.
    """
    assert isinstance(the_job, models.Job)
    logger.debug("In find dependent jobs for %s" % the_job)
    if growing_list is None:
        growing_list: List[models.Job] = []
    descendent_files = models.File.objects.filter(job=the_job)
    if not _follows_imported_files(the_job, follow_imported_files):
        descendent_files = descendent_files.exclude(fileimport__isnull=False)
    descendent_file_uses = models.FileUse.objects.filter(file__in=descendent_files)

    # Build set of dependent jobs, handling cases where jobs may have been deleted
    # (CASCADE delete removes FileUse when its Job is deleted)
    dependent_job_set = set()
    for file_use in descendent_file_uses:
        try:
            # Access the job - this may raise DoesNotExist if job was cascade-deleted
            if file_use.job and file_use.job != the_job:
                dependent_job_set.add(file_use.job)
        except models.Job.DoesNotExist:
            # Job was already deleted (cascade delete from parent), skip it
            logger.debug("Skipping FileUse with deleted job")
            continue

    unsorted_descendent_jobs = list(dependent_job_set) + list(models.Job.objects.filter(parent=the_job))
    # Reduce to uniques
    unsorted_descendent_jobs = list(set(unsorted_descendent_jobs))
    # Sort so that "leaves" will be processed first
    descendent_jobs = sorted(
        unsorted_descendent_jobs, key=version_sort_key, reverse=True
    )
    logger.debug(
        "descendent_jobs of %s: {%s}",
        the_job,
        [dj.number for dj in descendent_jobs],
    )
    original_length = len(growing_list)
    for descendent_job in descendent_jobs:
        if descendent_job not in growing_list:
            growing_list.append(descendent_job)
            find_dependent_jobs(
                descendent_job, growing_list, leaf_action, follow_imported_files
            )

    if original_length == len(growing_list):
        logger.debug("This node added no new jobs")
        # Here if (after performing leaf action on any relevant descendents) the list has grown
        if leaf_action is not None:
            leaf_action(the_job, growing_list)
    else:
        logger.debug("Growing list is now %s", [j.number for j in growing_list])

    return growing_list


def delete_job_and_dir(the_job: models.Job, growing_list: List[models.Job]):
    """Delete a job, every file it owns (imported ones included) and its directory."""
    logger.warning("Deleting job %s", the_job)
    for char_value_of_job in the_job.char_values.all():
        char_value_of_job.delete()
    for float_value_of_job in the_job.float_values.all():
        float_value_of_job.delete()
    job_file: models.File
    for job_file in the_job.files.all():
        try:
            job_file.path.unlink()
        except FileNotFoundError:
            logger.error("File  not found when trying to delete it %s", job_file.path)
        job_file.delete()
    logger.warning("Deleting directory %s", the_job.directory)
    if the_job.directory.exists() and the_job.directory.is_dir():
        shutil.rmtree(str(the_job.directory))
    logger.info("Deleted directory %s", the_job.directory)
    the_job.delete()
    if the_job in growing_list:
        growing_list.remove(the_job)


def _top_level_ancestor(job: models.Job) -> models.Job:
    while job.parent_id is not None:
        job = job.parent
    return job


def _empty_job_directory(directory: Path):
    if not (directory.exists() and directory.is_dir()):
        return
    for child in directory.iterdir():
        if child.name in FILE_HOLDER_KEEPS:
            continue
        if child.is_dir() and not child.is_symlink():
            shutil.rmtree(str(child))
        else:
            child.unlink()


def delete_job_keeping_imported_files(
    the_job: models.Job, growing_list: List[models.Job]
):
    """Delete ``the_job`` but keep the files it imported.

    A job that imported nothing is simply deleted. One that did stays as a
    FILE_HOLDER: its outputs, values, working files and its uses of anything
    else go, and its imported files (rows and bytes) stay, as do all uses of
    them. A
    sub-job's imported files are handed to its top-level job instead, which
    is either kept or becomes a file holder in its own turn: sub-jobs are
    deleted before their parent, which would otherwise cascade a sub-job
    file holder away.
    """
    kept_ids = list(imported_files_of(the_job).values_list("id", flat=True))
    if not kept_ids or the_job.status == models.Job.Status.FILE_HOLDER:
        delete_job_and_dir(the_job, growing_list)
        return

    if the_job.parent_id is not None:
        top = _top_level_ancestor(the_job)
        logger.warning(
            "Handing %d imported file(s) of job %s to job %s",
            len(kept_ids), the_job, top,
        )
        models.File.objects.filter(id__in=kept_ids).update(job=top)
        delete_job_and_dir(the_job, growing_list)
        return

    logger.warning(
        "Keeping job %s as a file holder for %d imported file(s)",
        the_job, len(kept_ids),
    )
    the_job.char_values.all().delete()
    the_job.float_values.all().delete()
    for job_file in the_job.files.exclude(id__in=kept_ids):
        try:
            job_file.path.unlink()
        except FileNotFoundError:
            logger.error("File  not found when trying to delete it %s", job_file.path)
        job_file.delete()
    # It no longer uses or makes anything but the files it holds; its report
    # lists those through its own uses of them. Uses of the kept files by
    # OTHER jobs are theirs, and stay.
    models.FileUse.objects.filter(job=the_job).exclude(
        file_id__in=kept_ids, role=models.FileUse.Role.IN
    ).delete()
    _empty_job_directory(the_job.directory)
    the_job.status = models.Job.Status.FILE_HOLDER
    the_job.evaluation = models.Job.Evaluation.UNKNOWN
    the_job.save(update_fields=["status", "evaluation"])
    if the_job in growing_list:
        growing_list.remove(the_job)


def deletion_leaf_action(
    delete_imported_files: bool, handled: Optional[Set[int]] = None
) -> Callable[[models.Job, List[models.Job]], None]:
    """The leaf action that deletes a job, by what happens to its imports.

    ``handled`` collects the id of every job acted on: a job kept as a file
    holder still exists afterwards, and must not be deleted a second time.
    """
    action = (
        delete_job_and_dir
        if delete_imported_files
        else delete_job_keeping_imported_files
    )
    if handled is None:
        return action

    def recording_action(job: models.Job, growing_list: List[models.Job]):
        handled.add(job.id)
        action(job, growing_list)

    return recording_action


def delete_job_and_dependents(
    the_job: models.Job, delete_imported_files: bool = False
):
    logger.warning("Deleting job %s and its dependents" % the_job)
    find_dependent_jobs(
        the_job,
        leaf_action=deletion_leaf_action(delete_imported_files),
        follow_imported_files=delete_imported_files,
    )


def find_bulk_dependent_jobs(
    job_ids: List[int], delete_imported_files: bool = False
) -> dict:
    """For a list of job IDs, find the union of all their dependent jobs,
    excluding jobs already in the selection.

    Returns a dict with:
      - selected_jobs: Job objects for the given IDs
      - additional_dependents: dependent jobs NOT in the selection
      - all_jobs_to_delete: union of selected + additional dependents
      - has_active_dependents: True if any additional dependent is running/queued
    """
    selected_jobs = list(models.Job.objects.filter(id__in=job_ids))
    selected_id_set = set(job_ids)

    all_dependents = set()
    for job in selected_jobs:
        deps = find_dependent_jobs(
            job, follow_imported_files=delete_imported_files
        )
        all_dependents.update(deps)

    additional_dependents = sorted(
        [j for j in all_dependents if j.id not in selected_id_set],
        key=version_sort_key,
    )

    has_active_dependents = any(
        j.status in ACTIVE_STATUSES for j in additional_dependents
    )

    return {
        "selected_jobs": selected_jobs,
        "additional_dependents": additional_dependents,
        "all_jobs_to_delete": selected_jobs + additional_dependents,
        "has_active_dependents": has_active_dependents,
    }


def imported_files_at_stake(
    jobs: Iterable[models.Job], selected_jobs: Iterable[models.Job]
) -> List[Dict]:
    """The imported files the delete-or-keep choice applies to, for ``jobs``
    (every job that goes if the imported files go), each with the jobs that
    used it --- other than the ``selected_jobs`` and their sub-jobs, which go
    either way. Those users go too if the files are deleted.

    A file holder's files are left out: deleting a file holder deletes them
    whichever way the choice goes.
    """
    jobs = list(jobs)
    job_ids = {j.id for j in jobs}
    selected_numbers = [j.number for j in selected_jobs]

    def is_selected(job: models.Job) -> bool:
        return any(
            job.number == n or job.number.startswith(n + ".")
            for n in selected_numbers
        )

    files = (
        models.File.objects.filter(job_id__in=job_ids, fileimport__isnull=False)
        .exclude(job__status=models.Job.Status.FILE_HOLDER)
        .select_related("job", "fileimport")
        .order_by("job_id", "id")
    )
    result = []
    for imported in files:
        users = [
            j
            for j in models.Job.objects.filter(
                file_uses__file=imported, file_uses__role=models.FileUse.Role.IN
            ).distinct()
            if not is_selected(j)
        ]
        result.append(
            {
                "id": imported.id,
                "uuid": str(imported.uuid),
                "name": imported.name,
                "source_name": Path(imported.fileimport.name).name,
                "annotation": imported.annotation,
                "job": _job_brief(imported.job),
                "used_by": [
                    _job_brief(j) for j in sorted(users, key=version_sort_key)
                ],
            }
        )
    return result


def _job_brief(job: models.Job) -> Dict:
    return {
        "id": job.id,
        "number": job.number,
        "title": job.title,
        "parent": job.parent_id,
        "status": job.status,
    }


def delete_multiple_jobs_and_dependents(
    job_ids: List[int], delete_imported_files: bool = False
):
    """Delete multiple jobs and all their dependents.

    Uses find_bulk_dependent_jobs to get the complete set, then deletes
    each job and its dependents leaf-first. Tracks already-handled jobs
    to avoid double-deletion when dependency graphs overlap (a job kept as a
    file holder is still there afterwards).
    """
    bulk_info = find_bulk_dependent_jobs(job_ids, delete_imported_files)
    all_to_delete = sorted(
        bulk_info["all_jobs_to_delete"], key=version_sort_key, reverse=True
    )

    handled: Set[int] = set()
    leaf_action = deletion_leaf_action(delete_imported_files, handled)
    for job in all_to_delete:
        if job.id in handled:
            continue
        try:
            job.refresh_from_db()
        except models.Job.DoesNotExist:
            handled.add(job.id)
            continue
        find_dependent_jobs(
            job,
            leaf_action=leaf_action,
            follow_imported_files=delete_imported_files,
        )
        handled.add(job.id)
