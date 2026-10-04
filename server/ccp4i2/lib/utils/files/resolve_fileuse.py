"""Resolve a fileUse reference to file metadata -- the contracted surface.

This is the adapter for ``GET /projects/{id}/resolve_fileuse/`` and for
``manage.py set_job_parameter``. The reference syntax and its resolution live in
:mod:`ccp4i2.lib.utils.files.file_use`, which the i2run command renderer and
i2run's own argument parser also use, so that the endpoint and the command line
cannot disagree about what ``[-1].XYZOUT`` means.

They did disagree, and worse than that: **nothing here worked at all.** Every
path built its return value with ``Result.failure(...)`` (nine places) or
``Result.success(...)`` (one), and ``Result`` implements neither -- it has
``ok`` and ``fail``. So the endpoint raised ``AttributeError`` on success and on
failure alike and could only ever 500, and the contract test passed because it
exercised the parser alone. The ordinal semantics the service contract
described were unreachable for the same reason, and would have been arbitrary
anyway: jobs were ordered by ``Job.number``, a CharField, so job 11 sorted
before job 2.

Two behaviours are preserved deliberately:

* the response shape, including ``fullPath``, which the published
  ``ResolveFileUseResponse`` type promises;
* a reference with no direction resolves against what a job PRODUCED first and
  what it consumed second, which is what the old code did by looking at File
  rows before FileUse rows. The i2run command line says the direction outright
  (``fileOut=`` / ``fileIn=``), which is better, but this endpoint's callers
  cannot, so it keeps the precedence.
"""

import logging

from ...response import Result
from .file_use import (
    FILE_IN,
    FILE_OUT,
    FileUseError,
    parse_file_use,
    resolve_file_reference,
)

logger = logging.getLogger(f"ccp4i2:{__name__}")


def is_fileuse_pattern(value: str) -> bool:
    """Does *value* look like a fileUse reference, with no ``fileUse=`` prefix?

    Used to auto-detect the syntax. Deliberately requires brackets: a bare
    ``3.XYZOUT`` is accepted by the parser as sugar on a command line, but
    auto-detecting it here would misread any ordinary dotted string.
    """
    if not value or not isinstance(value, str):
        return False
    if "[" not in value or "]" not in value or "." not in value:
        return False
    try:
        parse_file_use(value)
        return True
    except FileUseError:
        return False


def parse_fileuse(fileuse: str) -> dict:
    """Parse a fileUse string into the four contracted keys.

    ``{task_name, jobIndex, jobParamName, paramIndex}`` is the shape the
    service contract documents and consumers pattern-match on, so it is kept
    even though :func:`parse_file_use` carries more (an absolute job number is
    distinct from an ordinal there).

    Raises ValueError, as before; FileUseError is a ValueError.
    """
    if fileuse.startswith("fileUse="):
        fileuse = fileuse[len("fileUse=") :]
    for quote in ('"', "'"):
        if fileuse.startswith(quote) and fileuse.endswith(quote):
            fileuse = fileuse[1:-1]

    ref = parse_file_use(fileuse)
    return {
        "task_name": ref.task_name,
        # An absolute reference names a job number; an ordinal is an index.
        # Both surface here as jobIndex, which is what the contract has.
        "jobIndex": int(ref.job_number) if ref.job_number is not None else ref.index,
        "jobParamName": ref.param_name,
        "paramIndex": ref.param_index,
    }


def _as_uuid(text):
    import uuid
    try:
        return uuid.UUID(str(text).strip())
    except (ValueError, AttributeError):
        return None


def resolve_fileuse(project, fileuse: str):
    """Resolve *fileuse* against *project*, as ``Result.ok`` / ``Result.fail``.

    Produced files are consulted before consumed ones, so an unqualified
    reference keeps meaning what it did.
    """
    if fileuse.startswith("fileUse="):
        fileuse = fileuse[len("fileUse=") :]

    # A file's own id (as project_jobs and the file list give it) names it
    # unambiguously, including a file that went into a list element, which
    # a [job].PARAM reference cannot name.
    file_uuid = _as_uuid(fileuse)
    if file_uuid is not None:
        from ....db import models
        from .file_use import file_dict_for_file
        the_file = models.File.objects.filter(uuid=file_uuid, job__project=project).first()
        if the_file is None:
            return Result.fail(f"no file {fileuse} in this project")
        return Result.ok(file_dict_for_file(the_file, with_full_path=True))

    first_error = None
    for keyword in (FILE_OUT, FILE_IN):
        try:
            file_dict = resolve_file_reference(
                project, keyword, fileuse, with_full_path=True
            )
        except FileUseError as err:
            if first_error is None:
                first_error = err
            continue
        logger.debug("Resolved fileUse '%s' to %s", fileuse, file_dict["baseName"])
        return Result.ok(file_dict)

    return Result.fail(str(first_error))



def parse_fileuse_value(value: str) -> str:
    """Strip a ``fileUse=`` prefix and any surrounding quotes from *value*.

    ``manage.py set_job_parameter`` has imported this since it was written, and
    it never existed -- so that command failed at import with ``ImportError``,
    every invocation, not only the fileUse ones. Supplied here rather than
    dropped from the command, because the command wants it.
    """
    value = (value or "").strip()
    if value.startswith("fileUse="):
        value = value[len("fileUse=") :]
    for quote in ('"', "'"):
        if len(value) >= 2 and value.startswith(quote) and value.endswith(quote):
            value = value[1:-1]
    return value
