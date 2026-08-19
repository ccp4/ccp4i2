"""Move a project to a new location on disk, or repair its paths in place.

The CCP4i2 data model is very nearly location-independent: ``CDataFile`` is
serialised as ``project`` (UUID) + ``relPath`` + ``baseName``, and the only
absolute path in the database is :attr:`Project.directory`. Moving a project is
therefore mostly ``mv`` plus one field update.

What spoils it is *leakage*: a handful of places bake an absolute path into a
file inside the project directory.

* Programs that write their own XML (``program.xml`` from aimless/pointless
  carries ``<MTZmergedfilename>``). We cannot fix these at source.
* Generated scripts that are re-read by a viewer -- Coot's ``script.py`` and
  saved state, ccp4mg's ``*.scene.xml``.
* Legacy database dumps (``DATABASE.db.xml``, ``job.ccp4db.xml``) read on import.
* Cached reports, which are simply deleted here and regenerated on demand.

This module rewrites the first three groups and deletes the fourth. It
deliberately leaves log files, stdout/stderr and command scripts alone: those
are a historical record of what actually ran, and rewriting them would falsify
it.

Two entry points:

``move_project``
    Relocate the bytes and rebase the paths. Uses :func:`os.rename` when the
    destination is on the same filesystem (atomic, instant, no extra disk) and
    falls back to copy-verify-delete across devices.

``repair_project_paths``
    Rebase paths without moving anything, for the case where the directory
    stayed put but its absolute path changed underneath the user -- a renamed
    mount point, a reorganised network share.
"""

import errno
import json
import logging
import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from django.db import transaction

from .models import Job, Project

logger = logging.getLogger(f"ccp4i2:{__name__}")


#: Written at the destination root while a move is in flight so that a crashed
#: move can be diagnosed (and the half-rebased tree recognised) afterwards.
JOURNAL_NAME = ".ccp4i2-move.json"

#: Regenerated on demand by ``get_job_report_xml``; deleting is cheaper and
#: safer than rewriting.
CACHE_PATTERNS = (
    "report_xml.xml",
    "report.html",
    "report_tmp.html",
    "report_tmp.previous_*.html",
    "diagnostic_report.html",
)

#: Provenance. These record what ran, with the paths it ran against; rewriting
#: them would misrepresent history, and nothing reads them back as input.
PROVENANCE_NAMES = frozenset(
    {
        "log.txt",
        "stdout.txt",
        "stderr.txt",
        "com.txt",
        "diagnostic.xml",
    }
)
PROVENANCE_PREFIXES = ("stdout.", "stderr.", "log.")
PROVENANCE_SUFFIXES = (".log", ".txt")

#: Never opened as text, whatever their content. Coordinate and reflection
#: files are excluded on principle as well as on prudence: a path inside a PDB
#: REMARK is part of that file's own provenance, not a CCP4i2 reference.
BINARY_SUFFIXES = frozenset(
    {
        ".mtz", ".map", ".ccp4", ".mrc", ".pkl", ".pickle", ".npy", ".npz",
        ".png", ".jpg", ".jpeg", ".gif", ".bmp", ".tif", ".tiff", ".ico",
        ".pdf", ".zip", ".gz", ".bz2", ".xz", ".tar", ".tgz", ".7z",
        ".so", ".dylib", ".dll", ".exe", ".pyc", ".pyo",
        ".pdb", ".ent", ".cif", ".mmcif", ".sca", ".hkl", ".mmtf", ".bin",
    }
)

#: Files larger than this are not slurped into memory. Anything holding a
#: project path reference is orders of magnitude smaller.
MAX_REWRITE_BYTES = 64 * 1024 * 1024


class MoveProjectError(Exception):
    """Raised when a move or repair cannot be completed."""


@dataclass
class RewriteCandidate:
    """One file that contains the old root and would be rewritten."""

    path: Path
    occurrences: int

    def as_dict(self, relative_to: Path) -> Dict:
        try:
            name = str(self.path.relative_to(relative_to))
        except ValueError:
            name = str(self.path)
        return {"path": name, "occurrences": self.occurrences}


@dataclass
class RebasePlan:
    """What a rebase would do. Safe to compute without touching anything."""

    old_root: str
    new_root: str
    rewrites: List[RewriteCandidate] = field(default_factory=list)
    caches: List[Path] = field(default_factory=list)
    skipped_binary: List[Path] = field(default_factory=list)
    skipped_provenance: List[Path] = field(default_factory=list)
    skipped_too_large: List[Path] = field(default_factory=list)

    @property
    def total_occurrences(self) -> int:
        return sum(candidate.occurrences for candidate in self.rewrites)

    def as_dict(self, relative_to: Path) -> Dict:
        return {
            "old_root": self.old_root,
            "new_root": self.new_root,
            "rewrites": [c.as_dict(relative_to) for c in self.rewrites],
            "total_occurrences": self.total_occurrences,
            "caches_to_delete": [str(p.relative_to(relative_to)) for p in self.caches],
            "skipped": {
                "binary": len(self.skipped_binary),
                "provenance": len(self.skipped_provenance),
                "too_large": [str(p) for p in self.skipped_too_large],
            },
        }


# ---------------------------------------------------------------------------
# Path variants
# ---------------------------------------------------------------------------


def path_variants(root: str) -> Dict[str, str]:
    """Spellings of ``root`` that may appear inside project files, by form.

    A Windows root shows up both as ``C:\\Users\\me\\Proj`` and, because plenty
    of code normalises separators, as ``C:/Users/me/Proj``; inside a JSON string
    or a Python literal the backslashes may also be doubled.

    Keyed by form rather than returned as a list so that an old root and a new
    root pair up by *meaning*. Pairing positionally would mismatch the moment
    the two roots have different shapes -- moving a POSIX project onto a Windows
    share, say -- and silently corrupt whatever it touched.
    """
    root = root.rstrip("/\\")
    posix = root.replace("\\", "/")
    windows = root.replace("/", "\\")
    return {
        "native": root,
        "posix": posix,
        "windows": windows,
        "windows_escaped": windows.replace("\\", "\\\\"),
    }


def substitution_map(old_root: str, new_root: str) -> List[tuple]:
    """Ordered (old, new) pairs, longest old-spelling first.

    Longest first so that the doubled-backslash spelling is matched before the
    plain one can consume half of it.
    """
    old_variants = path_variants(old_root)
    new_variants = path_variants(new_root)

    pairs = []
    for form, old in old_variants.items():
        new = new_variants[form]
        if old and old != new and (old, new) not in pairs:
            pairs.append((old, new))

    pairs.sort(key=lambda pair: len(pair[0]), reverse=True)
    return pairs


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def is_provenance(path: Path) -> bool:
    """True for files that record what ran, and so must not be rewritten."""
    name = path.name
    if name in PROVENANCE_NAMES:
        return True
    if any(name.startswith(prefix) for prefix in PROVENANCE_PREFIXES):
        return True
    return path.suffix.lower() in PROVENANCE_SUFFIXES


def is_cache(path: Path) -> bool:
    """True for regenerable report caches, which are deleted rather than fixed."""
    from fnmatch import fnmatch

    return any(fnmatch(path.name, pattern) for pattern in CACHE_PATTERNS)


def read_text_if_textual(path: Path) -> Optional[str]:
    """Return the file's text, or ``None`` if it is binary or unreadable.

    Sniffs for a NUL byte rather than trusting the extension, so an unexpected
    binary with a ``.xml`` name cannot be corrupted.
    """
    try:
        with open(path, "rb") as handle:
            head = handle.read(8192)
            if b"\0" in head:
                return None
            handle.seek(0)
            raw = handle.read()
    except OSError as err:
        logger.warning("Could not read %s: %s", path, err)
        return None
    for encoding in ("utf-8", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return None


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------


def plan_rebase(root: Path, old_root: str, new_root: str) -> RebasePlan:
    """Walk ``root`` and work out what a rebase would change.

    Purely read-only, so it can back a dry-run in the UI.
    """
    plan = RebasePlan(old_root=old_root, new_root=new_root)
    substitutions = substitution_map(old_root, new_root)
    needles = [old for old, _ in substitutions]

    for path in sorted(root.rglob("*")):
        if path.is_dir() or path.is_symlink():
            continue
        if path.name == JOURNAL_NAME:
            continue
        if is_cache(path):
            plan.caches.append(path)
            continue
        if is_provenance(path):
            plan.skipped_provenance.append(path)
            continue
        if path.suffix.lower() in BINARY_SUFFIXES:
            plan.skipped_binary.append(path)
            continue
        try:
            size = path.stat().st_size
        except OSError:
            continue
        if size > MAX_REWRITE_BYTES:
            plan.skipped_too_large.append(path)
            continue

        text = read_text_if_textual(path)
        if text is None:
            plan.skipped_binary.append(path)
            continue
        occurrences = sum(text.count(needle) for needle in needles)
        if occurrences:
            plan.rewrites.append(RewriteCandidate(path=path, occurrences=occurrences))

    return plan


# ---------------------------------------------------------------------------
# Detecting roots left over from an earlier location
# ---------------------------------------------------------------------------

#: The top-level directories every CCP4i2 project has. An absolute path
#: containing one of them implies a project root: everything to its left.
PROJECT_MARKERS = ("CCP4_JOBS", "CCP4_IMPORTED_FILES", "CCP4_COOT", "CCP4_PROJECT_FILES")

STALE_ROOT_RE = re.compile(
    r"(?:[A-Za-z]:[\\/]|/)[^\s\"'<>|*?]*?"
    r"(?=[\\/]{1,2}(?:" + "|".join(PROJECT_MARKERS) + r")\b)"
)


def _normalise_root(text: str) -> str:
    """Tidy a root recovered from arbitrary text into a comparable form.

    Un-doubles escaped backslashes and strips the extra slashes a ``file:///``
    URL contributes, while leaving a UNC prefix (``\\\\server\\share``,
    ``//server/share``) intact -- those leading separators are part of the path.
    """
    is_unc = text.startswith("\\\\\\\\") or text.startswith("//")
    text = text.replace("\\\\", "\\")
    if not is_unc:
        text = re.sub(r"^/{2,}", "/", text)
    elif text.startswith("///"):
        # file:/// followed by a rooted POSIX path, not a UNC share.
        text = re.sub(r"^/{2,}", "/", text)
    return text.rstrip("/\\")


def find_stale_roots(root: Path, current_root: str) -> Dict[str, int]:
    """Count references to project roots other than ``current_root``.

    A project that has been moved or imported before may still carry paths from
    an even earlier location -- a previous move, a different machine, a mount
    point that has since been renamed. A rebase only knows about the one root
    it was asked to replace, so this reports the others rather than leaving the
    user to discover them when something silently fails to open.
    """
    current = _normalise_root(str(current_root))
    counts: Dict[str, int] = {}

    for path in sorted(root.rglob("*")):
        if path.is_dir() or path.is_symlink():
            continue
        if path.name == JOURNAL_NAME or is_cache(path) or is_provenance(path):
            continue
        if path.suffix.lower() in BINARY_SUFFIXES:
            continue
        try:
            if path.stat().st_size > MAX_REWRITE_BYTES:
                continue
        except OSError:
            continue
        text = read_text_if_textual(path)
        if text is None:
            continue
        for match in STALE_ROOT_RE.finditer(text):
            found = _normalise_root(match.group(0))
            if found and found != current:
                counts[found] = counts.get(found, 0) + 1

    return dict(sorted(counts.items(), key=lambda item: item[1], reverse=True))


# ---------------------------------------------------------------------------
# Applying
# ---------------------------------------------------------------------------


def apply_rebase(plan: RebasePlan, delete_caches: bool = True) -> List[Path]:
    """Carry out ``plan``. Returns the files actually rewritten.

    The caller keeps the returned list so that a later failure can be undone by
    running the inverse substitution over exactly those files.
    """
    substitutions = substitution_map(plan.old_root, plan.new_root)
    rewritten: List[Path] = []

    for candidate in plan.rewrites:
        text = read_text_if_textual(candidate.path)
        if text is None:
            continue
        updated = text
        for old, new in substitutions:
            updated = updated.replace(old, new)
        if updated == text:
            continue
        _write_preserving_mode(candidate.path, updated)
        rewritten.append(candidate.path)

    if delete_caches:
        for cache in plan.caches:
            try:
                cache.unlink()
            except OSError as err:
                logger.warning("Could not delete cached report %s: %s", cache, err)

    logger.info(
        "Rebased %d file(s) from %s to %s", len(rewritten), plan.old_root, plan.new_root
    )
    return rewritten


def undo_rebase(paths: List[Path], old_root: str, new_root: str) -> None:
    """Run the inverse substitution over ``paths``, best effort."""
    substitutions = substitution_map(new_root, old_root)
    for path in paths:
        text = read_text_if_textual(path)
        if text is None:
            continue
        updated = text
        for old, new in substitutions:
            updated = updated.replace(old, new)
        if updated != text:
            try:
                _write_preserving_mode(path, updated)
            except OSError as err:
                logger.error("Could not roll back rewrite of %s: %s", path, err)


def _write_preserving_mode(path: Path, text: str) -> None:
    """Write ``text`` to ``path``, keeping its permission bits."""
    try:
        mode = path.stat().st_mode
    except OSError:
        mode = None
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)
    if mode is not None:
        try:
            os.chmod(path, mode)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------


#: Statuses that mean a job may still be writing into the project directory.
#: PENDING is deliberately absent: it is the state of a task that has been
#: created and is still being edited, not one that has been submitted, so a
#: project full of half-filled-in tasks must still be movable. Matches the
#: active set used by ``api.views.active_jobs`` and
#: ``navigation.dependencies``.
ACTIVE_JOB_STATUSES = (
    Job.Status.QUEUED,
    Job.Status.RUNNING,
    Job.Status.RUNNING_REMOTELY,
)


def _running_jobs(project: Project) -> List[str]:
    return list(
        project.jobs.filter(status__in=ACTIVE_JOB_STATUSES).values_list(
            "number", flat=True
        )
    )


def preflight_move(project: Project, destination: Path) -> Path:
    """Validate a proposed move and return the resolved destination.

    Raises :class:`MoveProjectError` with a message fit to show the user.
    """
    source = Path(project.directory)
    destination = Path(destination).expanduser()

    if not source.is_dir():
        raise MoveProjectError(
            f"Project directory {source} does not exist. Use repair rather than move "
            "if the directory is intact but its path has changed."
        )

    # The server's working directory is not something the user can see, so a
    # relative destination would mean somewhere they did not choose.
    if not destination.is_absolute():
        raise MoveProjectError(
            f"Destination {destination} must be an absolute path."
        )
    destination = Path(os.path.normpath(str(destination)))

    if destination == source:
        raise MoveProjectError("Destination is the project's current location.")
    if destination.exists():
        raise MoveProjectError(f"Destination {destination} already exists.")

    parent = destination.parent
    if not parent.is_dir():
        raise MoveProjectError(f"Destination parent directory {parent} does not exist.")
    if not os.access(parent, os.W_OK):
        raise MoveProjectError(f"Destination parent directory {parent} is not writable.")

    resolved_source = source.resolve()
    if resolved_source in destination.resolve().parents or destination.resolve() == resolved_source:
        raise MoveProjectError(
            "Destination is inside the project directory; that would nest the "
            "project inside itself."
        )

    clash = (
        Project.objects.filter(directory=str(destination))
        .exclude(pk=project.pk)
        .first()
    )
    if clash is not None:
        raise MoveProjectError(
            f"Project '{clash.name}' is already registered at {destination}."
        )

    running = _running_jobs(project)
    if running:
        raise MoveProjectError(
            "Cannot move a project with jobs still running or queued "
            f"(job {', '.join(str(number) for number in running[:5])}). "
            "Wait for them to finish, or stop them first."
        )

    return destination


# ---------------------------------------------------------------------------
# Journal
# ---------------------------------------------------------------------------


def _write_journal(root: Path, old_root: str, new_root: str, state: str) -> None:
    try:
        (root / JOURNAL_NAME).write_text(
            json.dumps(
                {"old_root": old_root, "new_root": new_root, "state": state},
                indent=2,
            ),
            encoding="utf-8",
        )
    except OSError as err:
        logger.warning("Could not write move journal in %s: %s", root, err)


def _clear_journal(root: Path) -> None:
    try:
        (root / JOURNAL_NAME).unlink()
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------


def dry_run_move(project: Project, destination: Path) -> Dict:
    """Validate a move and report what it would rewrite, changing nothing."""
    destination = preflight_move(project, destination)
    source = Path(project.directory)
    plan = plan_rebase(source, str(source), str(destination))
    result = plan.as_dict(relative_to=source)
    result["source"] = str(source)
    result["destination"] = str(destination)
    result["same_filesystem"] = _same_filesystem(source, destination.parent)
    # Roots from an earlier life that this move will not touch.
    result["stale_roots"] = find_stale_roots(source, str(source))
    return result


def move_project(project: Project, destination: Path) -> Dict:
    """Move ``project`` to ``destination`` and rebase its paths.

    Same-filesystem moves use :func:`os.rename`: atomic, instant, and no
    second copy of what can be a multi-gigabyte project. Cross-device moves
    copy first and only delete the source once the database has been updated,
    so an interrupted move always leaves one complete, usable tree.
    """
    destination = preflight_move(project, destination)
    source = Path(project.directory)
    old_root, new_root = str(source), str(destination)

    cross_device = not _same_filesystem(source, destination.parent)
    logger.info(
        "Moving project '%s' from %s to %s (%s)",
        project.name,
        old_root,
        new_root,
        "copy across devices" if cross_device else "rename in place",
    )

    if cross_device:
        shutil.copytree(source, destination, symlinks=True)
    else:
        try:
            os.rename(source, destination)
        except OSError as err:
            # st_dev said same filesystem but rename disagreed -- bind mounts and
            # some network filesystems do this. Fall back rather than fail.
            if err.errno != errno.EXDEV:
                raise MoveProjectError(f"Could not move project: {err}") from err
            cross_device = True
            shutil.copytree(source, destination, symlinks=True)

    _write_journal(destination, old_root, new_root, "rebasing")

    rewritten: List[Path] = []
    try:
        plan = plan_rebase(destination, old_root, new_root)
        rewritten = apply_rebase(plan)

        with transaction.atomic():
            project.directory = new_root
            project.save(update_fields=["directory"])
    except Exception as err:
        logger.exception("Move failed after relocating bytes; rolling back")
        _rollback(source, destination, rewritten, old_root, new_root, cross_device)
        raise MoveProjectError(f"Move failed and was rolled back: {err}") from err

    _clear_journal(destination)

    if cross_device:
        try:
            shutil.rmtree(source)
        except OSError as err:
            # The move itself succeeded; the leftover is untidy, not broken.
            logger.error("Moved project but could not remove %s: %s", source, err)

    summary = plan.as_dict(relative_to=destination)
    summary.update(
        {
            "source": old_root,
            "destination": new_root,
            "rewritten": len(rewritten),
            "cross_device": cross_device,
            # Anything still pointing somewhere else entirely -- an earlier
            # move, an import from another machine, a renamed mount. Offer
            # these to the user as candidates for a follow-up repair.
            "stale_roots": find_stale_roots(destination, new_root),
        }
    )
    return summary


def repair_project_paths(
    project: Project, old_root: str, dry_run: bool = False
) -> Dict:
    """Rewrite stale ``old_root`` references without moving anything.

    For the case the directory never moved but its absolute path changed --
    a mount point renamed outside the user's control, a share reorganised.
    """
    new_root = str(Path(project.directory))
    old_root = str(old_root).rstrip("/\\")

    if not Path(new_root).is_dir():
        raise MoveProjectError(
            f"Project directory {new_root} does not exist, so there is nothing here "
            "to repair. If the directory is intact but its path changed, re-point "
            "the project at its current location instead."
        )
    if old_root == new_root:
        raise MoveProjectError(
            "The old and new roots are the same; there is nothing to repair."
        )

    plan = plan_rebase(Path(new_root), old_root, new_root)
    summary = plan.as_dict(relative_to=Path(new_root))
    summary.update({"source": old_root, "destination": new_root, "dry_run": dry_run})

    if dry_run:
        summary["rewritten"] = 0
        summary["stale_roots"] = find_stale_roots(Path(new_root), new_root)
        return summary

    rewritten = apply_rebase(plan)
    summary["rewritten"] = len(rewritten)
    summary["stale_roots"] = find_stale_roots(Path(new_root), new_root)
    return summary


# ---------------------------------------------------------------------------
# Re-pointing a project whose recorded directory has gone stale
# ---------------------------------------------------------------------------


def directory_is_missing(project: Project) -> bool:
    """True when the database points somewhere that is no longer there."""
    return not Path(project.directory).is_dir()


def looks_like_a_project(directory: Path) -> bool:
    """Cheap sanity check that ``directory`` really is a CCP4i2 project.

    Guards against a mistyped or mis-picked folder: re-pointing a project at
    somebody's Downloads folder and then rewriting paths through it would be a
    poor way to find out.
    """
    return any((directory / marker).is_dir() for marker in PROJECT_MARKERS)


def preflight_relocate(project: Project, new_directory: Path) -> Path:
    """Validate a re-point and return the resolved directory."""
    new_directory = Path(new_directory).expanduser()
    if not new_directory.is_absolute():
        raise MoveProjectError(f"{new_directory} must be an absolute path.")
    new_directory = Path(os.path.normpath(str(new_directory)))

    if not new_directory.is_dir():
        raise MoveProjectError(f"{new_directory} does not exist.")
    if not looks_like_a_project(new_directory):
        raise MoveProjectError(
            f"{new_directory} does not look like a CCP4i2 project directory: "
            "none of CCP4_JOBS, CCP4_IMPORTED_FILES, CCP4_COOT or "
            "CCP4_PROJECT_FILES is present."
        )
    if str(new_directory) == str(Path(project.directory)):
        raise MoveProjectError(
            "The project is already recorded at that location."
        )

    clash = (
        Project.objects.filter(directory=str(new_directory))
        .exclude(pk=project.pk)
        .first()
    )
    if clash is not None:
        raise MoveProjectError(
            f"Project '{clash.name}' is already registered at {new_directory}."
        )

    running = _running_jobs(project)
    if running:
        raise MoveProjectError(
            "Cannot re-point a project with jobs still running or queued "
            f"(job {', '.join(str(number) for number in running[:5])})."
        )

    return new_directory


def relocate_project(
    project: Project, new_directory: Path, dry_run: bool = False
) -> Dict:
    """Record that a project now lives somewhere else, and rebase its paths.

    Nothing is moved. This is for the case where the directory is intact but
    its absolute path changed underneath the user -- a drive renamed, a share
    remounted elsewhere, a projects folder moved with the Finder -- so the
    database is pointing at somewhere that no longer exists.

    The move operation cannot help there: it needs the files to still be where
    the database says they are.
    """
    new_directory = preflight_relocate(project, new_directory)
    old_root = str(Path(project.directory))
    new_root = str(new_directory)

    plan = plan_rebase(new_directory, old_root, new_root)
    summary = plan.as_dict(relative_to=new_directory)
    summary.update(
        {
            "source": old_root,
            "destination": new_root,
            "dry_run": dry_run,
            "directory_was_missing": not Path(old_root).is_dir(),
        }
    )

    if dry_run:
        summary["rewritten"] = 0
        summary["stale_roots"] = find_stale_roots(new_directory, new_root)
        return summary

    rewritten: List[Path] = []
    try:
        rewritten = apply_rebase(plan)
        with transaction.atomic():
            project.directory = new_root
            project.save(update_fields=["directory"])
    except Exception as err:
        logger.exception("Re-pointing project %s failed; rolling back", project.name)
        undo_rebase(rewritten, old_root, new_root)
        raise MoveProjectError(
            f"Could not re-point the project, and the rewrites were undone: {err}"
        ) from err

    logger.info("Re-pointed project '%s' from %s to %s", project.name, old_root, new_root)
    summary["rewritten"] = len(rewritten)
    summary["stale_roots"] = find_stale_roots(new_directory, new_root)
    return summary


# ---------------------------------------------------------------------------
# The whole store moved: one drive rename, every project stale
# ---------------------------------------------------------------------------


def _under(root: str, candidate: str) -> bool:
    """True when ``candidate`` is ``root`` itself or sits beneath it."""
    root = root.rstrip("/\\")
    if candidate == root:
        return True
    return candidate.startswith(root + "/") or candidate.startswith(root + "\\")


def plan_root_rebase(old_root: str, new_root: str) -> Dict:
    """Work out which projects a change of storage root accounts for.

    Renaming a drive or moving a projects folder invalidates every project at
    once, and re-pointing nineteen of them one at a time is no kind of recovery
    path. Given the old and new locations of the store, this maps each project
    beneath the old root onto its counterpart under the new one and reports
    which of them can actually be found there.
    """
    old_root = str(Path(old_root)).rstrip("/\\")
    new_root = str(Path(new_root)).rstrip("/\\")

    matched, missing, unaffected = [], [], []
    for project in Project.objects.all().order_by("name"):
        directory = str(Path(project.directory))
        if not _under(old_root, directory):
            unaffected.append({"id": project.pk, "name": project.name})
            continue
        remainder = directory[len(old_root):].lstrip("/\\")
        candidate = Path(new_root) / remainder if remainder else Path(new_root)
        entry = {
            "id": project.pk,
            "name": project.name,
            "old_directory": directory,
            "new_directory": str(candidate),
        }
        if candidate.is_dir() and looks_like_a_project(candidate):
            matched.append(entry)
        else:
            entry["reason"] = (
                "not found" if not candidate.is_dir() else "not a project directory"
            )
            missing.append(entry)

    return {
        "old_root": old_root,
        "new_root": new_root,
        "matched": matched,
        "missing": missing,
        "unaffected": unaffected,
    }


def rebase_projects_root(
    old_root: str, new_root: str, dry_run: bool = False
) -> Dict:
    """Re-point every project that moved with a change of storage root.

    Projects that cannot be found under the new root are reported and left
    alone; one project failing does not abandon the rest, since a partial
    recovery is better than none and the survivors can be retried.
    """
    summary = plan_root_rebase(old_root, new_root)
    summary["dry_run"] = dry_run

    if dry_run:
        summary["relocated"] = []
        summary["failed"] = []
        return summary

    relocated, failed = [], []
    for entry in summary["matched"]:
        project = Project.objects.get(pk=entry["id"])
        try:
            result = relocate_project(project, Path(entry["new_directory"]))
            relocated.append({**entry, "rewritten": result["rewritten"]})
        except Exception as err:
            logger.exception("Could not re-point project %s", project.name)
            failed.append({**entry, "error": str(err)})

    summary["relocated"] = relocated
    summary["failed"] = failed
    logger.info(
        "Re-pointed %d project(s) from %s to %s (%d failed, %d not found)",
        len(relocated),
        summary["old_root"],
        summary["new_root"],
        len(failed),
        len(summary["missing"]),
    )
    return summary


def update_projects_dir_preference(old_root: str, new_root: str) -> bool:
    """Follow the store move in ``preferences.json`` as well.

    Without this the database is repaired but the next new project is still
    created in the location that no longer exists.
    """
    from ..config import preferences

    prefs = preferences.load_preferences()
    current = prefs.get("projectsDir")
    if not current or not _under(str(Path(old_root)).rstrip("/\\"), str(Path(current))):
        return False

    remainder = str(Path(current))[len(str(Path(old_root)).rstrip("/\\")):].lstrip("/\\")
    updated = str(Path(new_root) / remainder) if remainder else str(Path(new_root))
    prefs["projectsDir"] = updated
    preferences.save_preferences(prefs)
    logger.info("Updated projectsDir preference from %s to %s", current, updated)
    return True


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _same_filesystem(source: Path, destination_parent: Path) -> bool:
    try:
        return os.stat(source).st_dev == os.stat(destination_parent).st_dev
    except OSError:
        return False


def _rollback(
    source: Path,
    destination: Path,
    rewritten: List[Path],
    old_root: str,
    new_root: str,
    cross_device: bool,
) -> None:
    """Put things back as they were, as far as we still can."""
    if cross_device:
        # The source was never touched -- discard the partial copy.
        shutil.rmtree(destination, ignore_errors=True)
        return

    undo_rebase(rewritten, old_root, new_root)
    _clear_journal(destination)
    try:
        os.rename(destination, source)
    except OSError as err:
        logger.error(
            "Could not rename %s back to %s: %s. The project directory is at %s; "
            "its database record still says %s.",
            destination,
            source,
            err,
            destination,
            source,
        )
