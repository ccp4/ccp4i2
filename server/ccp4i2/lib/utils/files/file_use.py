"""Naming a file on an i2run command line by where it came from.

A rendered command has to say which file each input is. A database id says it
exactly and says it unreadably: nobody knows, or can guess, that
``dbFileId=a3ed78ad466845a88765271c38d15149`` is the model from the last
refinement, and the first thing anyone does with a surfaced command is change
the inputs. A ``fileUse`` reference says the same thing in the terms people
actually think in:

    [3].XYZOUT                      XYZOUT of job 3
    [-1].XYZOUT                     XYZOUT of the most recent job
    prosmart_refmac[-1].XYZOUT      XYZOUT of the most recent prosmart_refmac
    prosmart_refmac[0].XYZOUT[1]    second XYZOUT of the first prosmart_refmac

**Brackets are not a matter of taste.** argparse rejects a bare negative:
``--XYZIN "-1.XYZOUT"`` fails with "expected at least one argument", because a
token starting with ``-`` is an option string and argparse's escape hatch only
exempts ``-1`` and ``-1.5``, not ``-1.XYZOUT``. So the bracketed form is the
only one in which relative indexing -- the whole point -- can be written.
A bare *positive* index (``3.XYZOUT``) is safe and accepted as sugar.

**What the integer means depends on whether a task is named**, and that is
deliberate:

* unqualified, ``>= 0``: a **job number**, as shown in the project. ``[3]`` is
  job 3, the thing a person reads off the screen.
* unqualified, ``< 0``: counting back from the most recent job.
* task-qualified: an **ordinal within that task's jobs**, so ``[0]`` is the
  first and ``[-1]`` the last. "Job number 2, but only if it happens to be a
  refmac" would be no use to anybody.

**Ordering is by creation, not by number.** ``Job.number`` is a CharField
holding a hierarchical number ("1", but also "1.1" for a sub-job), so sorting
it lexically puts job 11 before job 2, and there is no cast that orders "3.1"
sensibly among integers. Relative indexing therefore counts back through jobs
in creation order, which is what "the most recent" means anyway, and restricts
itself to top-level jobs -- the ones a person sees in the project list.

An absolute reference matches the job number exactly as a string, so a sub-job
can be named outright: ``[3.1].XYZOUT``.

**A relative reference skips jobs that have not got the file.** "The XYZOUT of
the most recent refmac" means the most recent refmac that produced one, not the
one someone configured a minute ago and has not run: a project accumulates
pending and failed jobs, and counting them would make ``[-1]`` resolve to
nothing useful most of the time. An *absolute* reference stays literal --- ask
for job 3 and you get job 3, or an error saying it has no such file.

The previous implementation indexed into an unordered queryset
(``jobs_query[jobIndex]``), so ``[-1]`` was not reliably the most recent --
which made relative indexing untrustworthy, and was the reason to write this
down.
"""

import difflib
import logging
import re
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(f"ccp4i2:{__name__}")


class FileUseError(ValueError):
    """A fileUse reference that cannot be parsed or cannot be resolved.

    Raised rather than returning None: the previous behaviour was to quietly
    drop the parameter, so ``--F_SIGF "fileUse=3.FREEROUT"`` configured an
    EMPTY plugin and the job failed later for an unrelated-looking reason.
    """


@dataclass(frozen=True)
class FileUseRef:
    """A parsed fileUse reference.

    Exactly one of *job_number* (an absolute reference, matched as a string so
    "3.1" works) and *index* (an ordinal, negative counting from the end) is
    set.

    *param_index* defaults to 0, the first file, as the service contract for
    ``GET /projects/{id}/resolve_fileuse/`` documents. The CLI and the endpoint
    must not disagree about what an omitted index means.
    """

    param_name: str
    task_name: Optional[str] = None
    job_number: Optional[str] = None
    index: Optional[int] = None
    param_index: int = 0
    #: The parameter text exactly as written, brackets included
    #: ("DICT_LIST[0]"). A real job_param_name can END in an index -- the File
    #: row for a CList element is recorded as 'DICT_LIST[0]' -- so "PARAM[n]"
    #: is genuinely ambiguous between a name and a name-plus-index. Resolution
    #: tries this literally before falling back to the split form.
    param_token: Optional[str] = None


# taskName[index].PARAM[index]
_QUALIFIED = re.compile(
    r"^(?P<task_name>[A-Za-z_]\w*)"
    r"\[(?P<index>-?\d+)\]"
    r"\.(?P<param_name>\w+)"
    r"(?:\[(?P<param_index>-?\d+)\])?$"
)
# [jobNumber].PARAM[index] or [-n].PARAM[index]
_BRACKETED = re.compile(
    r"^\[(?P<job_index>-?\d+(?:\.\d+)*)\]"
    r"\.(?P<param_name>\w+)"
    r"(?:\[(?P<param_index>-?\d+)\])?$"
)
# jobNumber.PARAM[index] -- positive only; see the module docstring
_BARE = re.compile(
    r"^(?P<job_index>\d+)"
    r"\.(?P<param_name>\w+)"
    r"(?:\[(?P<param_index>-?\d+)\])?$"
)


def _param_token(parts) -> str:
    """The parameter text as written: "DICT_LIST" or "DICT_LIST[0]"."""
    raw = parts.get("param_index")
    return (
        parts["param_name"] if raw is None else f"{parts['param_name']}[{raw}]"
    )


def parse_file_use(text: str) -> FileUseRef:
    """Parse a fileUse reference, or raise FileUseError saying what was wrong."""
    text = (text or "").strip()
    if not text:
        raise FileUseError("empty fileUse reference")

    def param_index_of(parts):
        raw = parts.get("param_index")
        return int(raw) if raw is not None else 0

    match = _QUALIFIED.match(text)
    if match:
        parts = match.groupdict()
        return FileUseRef(
            param_name=parts["param_name"],
            task_name=parts["task_name"],
            index=int(parts["index"]),
            param_index=param_index_of(parts),
            param_token=_param_token(parts),
        )

    for pattern in (_BRACKETED, _BARE):
        match = pattern.match(text)
        if not match:
            continue
        parts = match.groupdict()
        raw = parts["job_index"]
        if raw.startswith("-"):
            return FileUseRef(
                param_name=parts["param_name"],
                index=int(raw),
                param_index=param_index_of(parts),
                param_token=_param_token(parts),
            )
        return FileUseRef(
            param_name=parts["param_name"],
            job_number=raw,
            param_index=param_index_of(parts),
            param_token=_param_token(parts),
        )

    # A bare negative parses as a reference but can never reach us through
    # argparse, so say that rather than "bad syntax".
    if re.match(r"^-\d+\.\w+", text):
        index, _, rest = text.partition(".")
        raise FileUseError(
            f"'{text}': a negative index must be bracketed, as "
            f"'[{index}].{rest}' -- argparse treats a bare leading '-' as an "
            f"option, so this form can never be passed on a command line"
        )

    raise FileUseError(
        f"'{text}' is not a fileUse reference. Expected [jobNumber].PARAM, "
        f"[-n].PARAM, taskName[-n].PARAM or jobNumber.PARAM"
    )


def _did_you_mean(name: str, candidates) -> str:
    """" Did you mean 'X'?" for the nearest candidate, else "".

    A typo is the likeliest reason a reference does not resolve, and the
    alternative -- printing i2run's usage -- is no help at all here: a task like
    servalcat_pipe has 219 arguments, so the answer would be buried. The
    relevant list is always short (the parameters of one job, or the registered
    task names), so name it.
    """
    candidates = [str(c) for c in candidates if c]
    if not candidates:
        return ""
    close = difflib.get_close_matches(name, candidates, n=1, cutoff=0.6)
    if close and close[0] != name:
        return f" Did you mean '{close[0]}'?"
    # difflib is strict about case and about short strings; fall back to a
    # case-insensitive exact hit, which is a typo people make constantly.
    lowered = {c.lower(): c for c in candidates}
    if name.lower() in lowered and lowered[name.lower()] != name:
        return f" Did you mean '{lowered[name.lower()]}'?"
    return ""


def _and_these_exist(label: str, names, limit: int = 12) -> str:
    """A short, sorted inventory of what is actually available."""
    unique = sorted({str(n) for n in names if n})
    if not unique:
        return ""
    shown = unique[:limit]
    more = "" if len(unique) == len(shown) else f", and {len(unique) - len(shown)} more"
    return f" {label}: {', '.join(shown)}{more}."


def _known_task_names():
    from ....core.tasks import TASKS

    return list(TASKS)


def _param_names_on(job, role=None):
    """Parameter names *job* has a file under, optionally for one role only."""
    from ....db import models

    uses = models.FileUse.objects.filter(job=job)
    if role is not None:
        uses = uses.filter(role=role)
    names = set(uses.values_list("job_param_name", flat=True))
    if not names and role is None:
        names = set(
            models.File.objects.filter(job=job).values_list(
                "job_param_name", flat=True
            )
        )
    return names


def _candidate_jobs(project, ref: FileUseRef):
    """The jobs a reference could mean, in creation order.

    Creation order, not number order: see the module docstring.
    """
    from ....db import models

    jobs = models.Job.objects.filter(project=project)
    if ref.task_name is not None:
        known = _known_task_names()
        if ref.task_name not in known:
            raise FileUseError(
                f"'{ref.task_name}' is not a task."
                f"{_did_you_mean(ref.task_name, known)}"
            )
        return list(jobs.filter(task_name=ref.task_name).order_by("id"))
    # Unqualified: sub-jobs are not what anyone means by "the last job".
    return [job for job in jobs.order_by("id") if "." not in job.number]


def _role_of_owned(the_file, job):
    """Which direction *job* saw a file it owns, from its FileUse row."""
    from ....db import models

    use = models.FileUse.objects.filter(file=the_file, job=job).first()
    if use is not None:
        return use.role
    # No row: an imported file is an input to the job it was imported for,
    # anything else a job owns, it produced.
    return (
        models.FileUse.Role.IN
        if the_file.directory == models.File.Directory.IMPORT_DIR
        else models.FileUse.Role.OUT
    )


def _files_on(job, ref: FileUseRef, role):
    """Files *job* names as the reference's parameter, in the given direction.

    Two record types have to be consulted, and neither alone is enough:

    * **File** rows are the files *job* owns -- produced or imported. They carry
      the full parameter name, including a list index: a CList element is
      recorded as ``DICT_LIST[0]``.
    * **FileUse** rows are what *job* consumed, including files another job
      owns. ``fileIn=[2].FREERFLAG`` needs this: the file belongs to job 1.

    And the parameter text is tried literally ("DICT_LIST[0]") before the split
    form ("DICT_LIST", index 0), because a genuine parameter name can end in an
    index and the syntax cannot tell the two apart by itself.

    Beware a live-data wrinkle: for a CList element the File row says
    ``DICT_LIST[0]`` while the FileUse row for the very same file says only
    ``[0]`` -- the list's name is dropped when the use is recorded. That is a
    bug in the recording, not here, and it cannot be worked around safely by
    name (two file lists on one job would both record ``[0]``). Consulting File
    rows first is what makes list elements resolvable at all.

    Matched case-insensitively: no two parameters anywhere in the registry
    differ only by case. Never upper-cased -- capitalisation is a convention
    with real exceptions (PHIL parameters are lower-case, some classic ones are
    mixed, every task has jobTitle/jobStatus).
    """
    from ....db import models

    tokens = [ref.param_token or ref.param_name]
    if ref.param_name not in tokens:
        tokens.append(ref.param_name)

    for token in tokens:
        for lookup in (
            {"job_param_name": token},
            {"job_param_name__iexact": token},
        ):
            # Files this job owns, filtered to the direction it saw them.
            owned = [
                f
                for f in models.File.objects.filter(job=job, **lookup)
                .select_related("job__project")
                .order_by("id")
                if _role_of_owned(f, job) == role
            ]
            if owned:
                return owned

            # Files this job consumed, whoever owns them.
            file_ids = models.FileUse.objects.filter(
                job=job, role=role, **lookup
            ).values_list("file_id", flat=True)
            used = list(
                models.File.objects.filter(id__in=file_ids)
                .select_related("job__project")
                .order_by("id")
            )
            if used:
                return used
    return []


def _pick(files, ref: FileUseRef, what: str, job=None, role=None):
    if not files:
        available = _param_names_on(job, role) if job is not None else ()
        detail = (
            _did_you_mean(ref.param_name, available)
            + _and_these_exist("It has", available)
            if available
            else " That job has no files at all - has it run?"
        )
        raise FileUseError(f"{what} has no file for '{ref.param_name}'.{detail}")
    index = 0 if len(files) == 1 else ref.param_index
    try:
        return files[index]
    except IndexError:
        raise FileUseError(
            f"{what} has {len(files)} file(s) for '{ref.param_name}'; "
            f"index {ref.param_index} is out of range"
        ) from None


def _resolve_ref(project, ref: FileUseRef, role):
    """The file *ref* names in *project*."""
    from ....db import models

    if ref.job_number is not None:
        job = (
            models.Job.objects.filter(project=project, number=ref.job_number)
            .first()
        )
        if job is None:
            raise FileUseError(
                f"no job numbered {ref.job_number} in project '{project.name}'"
            )
        return _pick(
            _files_on(job, ref, role),
            ref,
            f"job {job.number} ({job.task_name})",
            job,
            role,
        )

    candidates = _candidate_jobs(project, ref)
    described = (
        f"{ref.task_name} jobs" if ref.task_name else "top-level jobs"
    )
    if not candidates:
        raise FileUseError(f"project '{project.name}' has no {described}")

    if ref.index is not None and ref.index < 0:
        # Relative: only jobs that actually have the file count.
        with_file = [
            (job, files)
            for job, files in (
                (job, _files_on(job, ref, role)) for job in candidates
            )
            if files
        ]
        if not with_file:
            available = set()
            for candidate in candidates:
                available |= _param_names_on(candidate, role)
            # An empty inventory is a different diagnosis from a wrong name:
            # the jobs exist but have produced nothing, which usually means
            # they have not been run.
            detail = (
                _did_you_mean(ref.param_name, available)
                + _and_these_exist("Between them they have", available)
                if available
                else f" None of those {len(candidates)} job(s) have any files "
                f"-- have they run?"
            )
            raise FileUseError(
                f"no {described} in project '{project.name}' have a file for "
                f"'{ref.param_name}'.{detail}"
            )
        try:
            job, files = with_file[ref.index]
        except IndexError:
            raise FileUseError(
                f"[{ref.index}]: only {len(with_file)} {described} in "
                f"'{project.name}' have a '{ref.param_name}'"
            ) from None
        return _pick(files, ref, f"job {job.number} ({job.task_name})", job, role)

    # Qualified, non-negative: a literal ordinal within that task's jobs.
    try:
        job = candidates[ref.index]
    except IndexError:
        raise FileUseError(
            f"{ref.task_name}[{ref.index}]: project '{project.name}' has "
            f"{len(candidates)} {described}"
        ) from None
    return _pick(
        _files_on(job, ref, role),
        ref,
        f"job {job.number} ({job.task_name})",
        job,
        role,
    )


def file_dict_for_file(the_file, with_full_path: bool = False) -> dict:
    """The key=value fields that identify *the_file* to a CDataFile.

    ``relPath`` comes from the job's own directory relative to the project's,
    not from a constructed ``CCP4_JOBS/job_<number>``: a sub-job numbered "3.1"
    lives at ``CCP4_JOBS/job_3/job_1``, which the constructed form gets wrong.

    ``fullPath`` is included only on request, because it must NOT reach a
    CDataFile: setting it alongside baseName rewrites baseName to the whole
    absolute path, which is how ``<baseName>/Users/.../xyzout.pdb</baseName>``
    appeared in a configured job. The service contract for
    ``GET /projects/{id}/resolve_fileuse/`` promises the field, so the endpoint
    asks for it; the i2run populator does not.
    """
    from pathlib import Path

    from ....db import models

    project_dir = Path(the_file.job.project.directory)
    file_dict = {
        "project": str(the_file.job.project.uuid).replace("-", ""),
        "baseName": the_file.name,
        "dbFileId": str(the_file.uuid).replace("-", ""),
    }
    # What the file holds, as the app sets it when a project file is picked
    # (set_input_by_context._set_file_from_db). Without it a file named by
    # fileIn=/fileOut= reached a task with contentFlag unset: a pipeline took
    # it, and the sub-job it handed it to refused it ("got 0, requires one of
    # IPAIR, FPAIR, IMEAN, FMEAN"), as phaser_simple_phil did.
    if the_file.content is not None:
        file_dict["contentFlag"] = the_file.content
    if the_file.sub_type is not None:
        file_dict["subType"] = the_file.sub_type
    if the_file.annotation:
        file_dict["annotation"] = the_file.annotation

    if the_file.directory == models.File.Directory.IMPORT_DIR:
        file_dict["relPath"] = "CCP4_IMPORTED_FILES"
        full_path = project_dir / "CCP4_IMPORTED_FILES" / the_file.name
    else:
        job_dir = Path(the_file.job.directory)
        try:
            file_dict["relPath"] = str(job_dir.relative_to(project_dir))
        except ValueError:
            # Not under the project (a relocated or imported project): keep
            # whatever of the path is meaningful rather than inventing one.
            parts = job_dir.parts
            if "CCP4_JOBS" in parts:
                file_dict["relPath"] = str(
                    Path(*parts[parts.index("CCP4_JOBS"):])
                )
            else:
                file_dict["relPath"] = job_dir.name
        full_path = job_dir / the_file.name

    if with_full_path:
        file_dict["fullPath"] = str(full_path)
    return file_dict


FILE_IN = "fileIn"
FILE_OUT = "fileOut"

#: Accepted but deprecated. ``fileUse=`` is what the CLI README documented and
#: what Qt-era i2run took, so scripts in the wild use it. It names no direction,
#: so it resolves against what a job produced before what it consumed -- the
#: same precedence the resolve_fileuse endpoint keeps, for the same reason.
FILE_USE = "fileUse"

FILE_KEYWORDS = (FILE_IN, FILE_OUT, FILE_USE)


def _role_of(keyword: str):
    """The FileUse role a keyword names, or None for the directionless alias."""
    from ....db import models

    if keyword == FILE_IN:
        return models.FileUse.Role.IN
    if keyword == FILE_OUT:
        return models.FileUse.Role.OUT
    if keyword == FILE_USE:
        return None
    raise FileUseError(
        f"'{keyword}=' is not a file reference. Use '{FILE_OUT}=' for a file a "
        f"job produced, or '{FILE_IN}=' for one it consumed"
    )


def resolve_db_file_id(file_id: str) -> dict:
    """The fields that identify the File with database id *file_id* to a
    CDataFile, as :func:`file_dict_for_file` gives them for a fileIn= /
    fileOut= reference. A bare dbFileId= on the command line needs them too:
    alone it left the file with no path when the job was validated, so its
    requiredContentFlag was never checked."""
    from ....db import models

    try:
        the_file = models.File.objects.get(uuid=file_id)
    except (models.File.DoesNotExist, ValueError):
        raise FileUseError(f"dbFileId={file_id}: no such file") from None
    return file_dict_for_file(the_file)


def resolve_file_reference(
    project, keyword: str, text: str, with_full_path: bool = False
) -> dict:
    """Resolve ``fileIn=``/``fileOut=`` *text* against *project*.

    *project* is a Project instance, its uuid, or its name. *with_full_path*
    adds ``fullPath``, which the resolve_fileuse endpoint's contract promises
    and which must not be passed to a CDataFile (see
    :func:`file_dict_for_file`).
    """
    from ....db import models

    if isinstance(project, models.Project):
        the_project = project
    else:
        text_project = str(project)
        try:
            the_project = models.Project.objects.get(uuid=text_project)
        except (models.Project.DoesNotExist, ValueError, TypeError):
            try:
                the_project = models.Project.objects.get(name=text_project)
            except models.Project.DoesNotExist:
                raise FileUseError(f"no project '{text_project}'") from None

    from ....db import models

    role = _role_of(keyword)
    ref = parse_file_use(text)

    if role is None:
        # Directionless (the fileUse= alias): produced first, then consumed.
        # Keeping one precedence rule here means the endpoint and the command
        # line cannot disagree about what an undirected reference means.
        first_error = None
        for candidate in (models.FileUse.Role.OUT, models.FileUse.Role.IN):
            try:
                the_file = _resolve_ref(the_project, ref, candidate)
            except FileUseError as err:
                if first_error is None:
                    first_error = err
                continue
            break
        else:
            raise first_error
    else:
        the_file = _resolve_ref(the_project, ref, role)

    logger.debug("%s=%s -> %s", keyword, text, the_file.name)
    return file_dict_for_file(the_file, with_full_path=with_full_path)


def file_reference_for_file(the_file, rendering_job=None):
    """``(keyword, reference)`` naming *the_file*, or None if there isn't one.

    A file is named by the job that OWNS it -- the one that produced it, or
    imported it -- which is what its File row records, and in the DIRECTION
    that job saw it, which is what its FileUse role records. So the MTZ
    imported for job 2 renders as ``fileIn=[2].F_SIGF`` to every later job that
    uses it, and job 2's own FREEROUT renders as ``fileOut=[2].FREEROUT``.

    The one case with no reference is the owning job itself: rendering job 2's
    command, ``[2].F_SIGF`` would point at the job being described, which is
    circular and useless to edit. Pass *rendering_job* to get None there, so
    the caller can fall back to the database id.

    (Testing the import directory instead, as this first did, was too blunt: it
    refused the useful case as well, and imported files are exactly the ones a
    person most wants to redirect at new data.)
    """
    from ....db import models

    if the_file.job is None or not the_file.job_param_name:
        return None
    if rendering_job is not None and the_file.job_id == rendering_job.id:
        return None

    use = models.FileUse.objects.filter(file=the_file, job=the_file.job).first()
    if use is not None:
        keyword = FILE_IN if use.role == models.FileUse.Role.IN else FILE_OUT
    else:
        # No row for the owning job: an import is an input to it, anything
        # else it owns, it produced.
        keyword = (
            FILE_IN
            if the_file.directory == models.File.Directory.IMPORT_DIR
            else FILE_OUT
        )
    return keyword, f"[{the_file.job.number}].{the_file.job_param_name}"
