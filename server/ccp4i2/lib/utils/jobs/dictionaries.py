"""Which ligand dictionaries belong with a job.

A job that is going to be drawn (a refinement, SubstituteLigand, a model
building session) has its ligand dictionary as one of its inputs or as one
of its outputs. That is the only association that is a fact rather than a
guess: a project routinely holds several ligands with the same residue name
(LIG, DRG), so "any dictionary in the project that defines LIG" is wrong as
soon as there are two.

The database already records both halves:

* a dictionary the job produced or imported is a ``File`` whose ``job`` is
  this job;
* a dictionary the job took from another job is a ``FileUse`` with role IN.

A subjob of a pipeline inherits its pipeline's dictionaries when it has none
of its own, because the pipeline is the job the user ran.
"""

from ccp4i2.db import models

DICTIONARY_TYPE = "application/refmac-dictionary"


def _own(job):
    return models.File.objects.filter(job=job, type__name=DICTIONARY_TYPE).order_by("id")


def _inputs(job):
    return (
        models.File.objects.filter(
            file_uses__job=job,
            file_uses__role=models.FileUse.Role.IN,
            type__name=DICTIONARY_TYPE,
        )
        .distinct()
        .order_by("id")
    )


def job_dictionaries(job):
    """``[(File, role)]`` for ``job``: its own dictionary files (role
    ``"own"``) then the ones it used as input (role ``"input"``), without
    duplicates. Falls back to the nearest ancestor that has any."""
    current = job
    while current is not None:
        seen = set()
        found = []
        for the_file in _own(current):
            seen.add(the_file.id)
            found.append((the_file, "own"))
        for the_file in _inputs(current):
            if the_file.id not in seen:
                seen.add(the_file.id)
                found.append((the_file, "input"))
        if found:
            return found
        current = current.parent
    return []


def companion_dictionaries(the_file):
    """The dictionaries that belong with a coordinate (or any) file: those of
    the job the file belongs to."""
    if the_file.job_id is None:
        return []
    return job_dictionaries(the_file.job)
