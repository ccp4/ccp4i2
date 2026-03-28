# Copyright (C) 2025-2026 University of York
# Copyright (C) 2025-2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
from typing import List
import logging
import shutil

from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")
logger.setLevel(logging.WARNING)


# Note that python handles tuple sort exactly as one would wish
# https://stackoverflow.com/questions/2574080/how-to-sort-in-python-a-list-of-strings-containing-numbers
def version_sort_key(job: models.Job) -> tuple:
    return tuple(map(int, job.number.split(".")))


def find_dependent_jobs(
    the_job: models.Job, growing_list: List[models.Job] = None, leaf_action=None
) -> List[models.Job]:
    assert isinstance(the_job, models.Job)
    logger.debug("In find dependent jobs for %s" % the_job)
    if growing_list is None:
        growing_list: List[models.Job] = []
    descendent_files = models.File.objects.filter(job=the_job)
    # Do I also have to worry about FileImports in this job ? I don't think so...those should have an associated FileUse
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
            find_dependent_jobs(descendent_job, growing_list, leaf_action)

    if original_length == len(growing_list):
        logger.debug("This node added no new jobs")
        # Here if (after performing leaf action on any relevant descendents) the list has grown
        if leaf_action is not None:
            leaf_action(the_job, growing_list)
    else:
        logger.debug("Growing list is now %s", [j.number for j in growing_list])

    return growing_list


def delete_job_and_dir(the_job: models.Job, growing_list: List[models.Job]):
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


def delete_job_and_dependents(the_job: models.Job):
    logger.warning("Deleting job %s and its dependents" % the_job)
    jobs_before = models.Job.objects.count()
    files_before = models.File.objects.count()
    file_uses_before = models.FileUse.objects.count()
    file_imports_before = models.FileImport.objects.count()
    find_dependent_jobs(the_job, leaf_action=delete_job_and_dir)
    jobs_after = models.Job.objects.count()
    files_after = models.File.objects.count()
    file_uses_after = models.FileUse.objects.count()
    file_imports_after = models.FileImport.objects.count()
    logger.info("Deleted %s jobs", jobs_after - jobs_before)
    logger.info("Deleted %s files", files_after - files_before)
    logger.info("Deleted %s file_uses", file_uses_after - file_uses_before)
    logger.info("Deleted %s jobs", file_imports_after - file_imports_before)


def find_bulk_dependent_jobs(job_ids: List[int]) -> dict:
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
        deps = find_dependent_jobs(job)
        all_dependents.update(deps)

    additional_dependents = sorted(
        [j for j in all_dependents if j.id not in selected_id_set],
        key=version_sort_key,
    )

    ACTIVE_STATUSES = {
        models.Job.Status.QUEUED,
        models.Job.Status.RUNNING,
        models.Job.Status.RUNNING_REMOTELY,
    }
    has_active_dependents = any(
        j.status in ACTIVE_STATUSES for j in additional_dependents
    )

    return {
        "selected_jobs": selected_jobs,
        "additional_dependents": additional_dependents,
        "all_jobs_to_delete": selected_jobs + additional_dependents,
        "has_active_dependents": has_active_dependents,
    }


def delete_multiple_jobs_and_dependents(job_ids: List[int]):
    """Delete multiple jobs and all their dependents.

    Uses find_bulk_dependent_jobs to get the complete set, then deletes
    each job and its dependents leaf-first. Tracks already-deleted jobs
    to avoid double-deletion when dependency graphs overlap.
    """
    bulk_info = find_bulk_dependent_jobs(job_ids)
    all_to_delete = sorted(
        bulk_info["all_jobs_to_delete"], key=version_sort_key, reverse=True
    )

    deleted_ids = set()
    for job in all_to_delete:
        if job.id not in deleted_ids:
            try:
                job.refresh_from_db()
            except models.Job.DoesNotExist:
                deleted_ids.add(job.id)
                continue
            find_dependent_jobs(job, leaf_action=delete_job_and_dir)
            deleted_ids.add(job.id)
