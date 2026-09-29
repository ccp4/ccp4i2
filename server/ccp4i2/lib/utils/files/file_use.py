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
    """

    param_name: str
    task_name: Optional[str] = None
    job_number: Optional[str] = None
    index: Optional[int] = None
    param_index: int = -1


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


def parse_file_use(text: str) -> FileUseRef:
    """Parse a fileUse reference, or raise FileUseError saying what was wrong."""
    text = (text or "").strip()
    if not text:
        raise FileUseError("empty fileUse reference")

    def param_index_of(parts):
        raw = parts.get("param_index")
        return int(raw) if raw is not None else -1

    match = _QUALIFIED.match(text)
    if match:
        parts = match.groupdict()
        return FileUseRef(
            param_name=parts["param_name"],
            task_name=parts["task_name"],
            index=int(parts["index"]),
            param_index=param_index_of(parts),
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
            )
        return FileUseRef(
            param_name=parts["param_name"],
            job_number=raw,
            param_index=param_index_of(parts),
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


def _files_on(job, param_name: str, role):
    """Files *job* names as *param_name* in the given direction.

    The role is the whole point of splitting the syntax into ``fileIn=`` and
    ``fileOut=``. Four tasks carry the same file parameter name in both
    inputData and outputData -- ``molrep_pipe`` and
    ``dr_mr_modelbuild_pipeline`` pass ``F_SIGF`` and ``FREERFLAG`` straight
    through, ``servalcat_pipe`` has ``METALCOORD_RESTRAINTS``,
    ``adding_stats_to_mmcif_i2`` has ``FPHIOUT``/``DIFFPHIOUT`` -- so a single
    keyword would have needed an arbitrary precedence rule to decide what
    ``[5].F_SIGF`` meant. Naming the direction removes the ambiguity instead of
    resolving it by fiat.

    Matched case-insensitively, because getting the case wrong is the typo
    people actually make and there is nothing to gain by refusing it: no two
    parameters anywhere in the registry differ only by case (7496 upper, 2036
    lower, 249 mixed, zero collisions). What this must NOT do is upper-case the
    name: parameter names are only *conventionally* capitalised, and
    PHIL-derived ones are lower-case
    (``phaser__crystal_symmetry__unit_cell``), a few classic ones are mixed
    (``F_SIGFanom``), and every task carries ``jobTitle``/``jobStatus``.
    """
    from ....db import models

    for lookup in (
        {"job_param_name": param_name},
        {"job_param_name__iexact": param_name},
    ):
        file_ids = models.FileUse.objects.filter(
            job=job, role=role, **lookup
        ).values_list("file_id", flat=True)
        files = list(
            models.File.objects.filter(id__in=file_ids)
            .select_related("job__project")
            .order_by("id")
        )
        if files:
            return files

        # Defensive fallback for a file with no FileUse row for its own job.
        # An imported file is an INPUT to the job it was imported for; anything
        # else a job owns, it produced.
        want_import = role == models.FileUse.Role.IN
        owned = [
            f
            for f in models.File.objects.filter(job=job, **lookup)
            .select_related("job__project")
            .order_by("id")
            if (f.directory == models.File.Directory.IMPORT_DIR) == want_import
        ]
        if owned:
            return owned
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
    try:
        return files[ref.param_index]
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
            _files_on(job, ref.param_name, role),
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
                (job, _files_on(job, ref.param_name, role)) for job in candidates
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
        _files_on(job, ref.param_name, role),
        ref,
        f"job {job.number} ({job.task_name})",
        job,
        role,
    )


def file_dict_for_file(the_file) -> dict:
    """The key=value fields that identify *the_file* to a CDataFile."""
    # Imported here, not at module scope: parse_file_use is pure, and a
    # pure parser should be testable without a configured Django.
    from ....db import models

    file_dict = {
        "project": str(the_file.job.project.uuid).replace("-", ""),
        "baseName": the_file.name,
        "dbFileId": str(the_file.uuid).replace("-", ""),
    }
    if the_file.directory == models.File.Directory.IMPORT_DIR:
        file_dict["relPath"] = "CCP4_IMPORTED_FILES"
    else:
        file_dict["relPath"] = f"CCP4_JOBS/job_{the_file.job.number}"
    return file_dict


FILE_IN = "fileIn"
FILE_OUT = "fileOut"


def _role_of(keyword: str):
    """The FileUse role a keyword names."""
    from ....db import models

    if keyword == FILE_IN:
        return models.FileUse.Role.IN
    if keyword == FILE_OUT:
        return models.FileUse.Role.OUT
    raise FileUseError(
        f"'{keyword}=' is not a file reference. Use '{FILE_OUT}=' for a file a "
        f"job produced, or '{FILE_IN}=' for one it consumed"
    )


def resolve_file_reference(project, keyword: str, text: str) -> dict:
    """Resolve ``fileIn=``/``fileOut=`` *text* against *project*.

    *project* is a Project instance, its uuid, or its name.
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

    role = _role_of(keyword)
    the_file = _resolve_ref(the_project, parse_file_use(text), role)
    logger.info("%s=%s -> %s", keyword, text, the_file.name)
    return file_dict_for_file(the_file)


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
