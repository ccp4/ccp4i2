"""
Structural pre-flight analysis for legacy CCP4i2 migration.

Given the set of legacy projects (id, name, on-disk directory) plus the migration
mode (schlurp vs copy) and a destination root, this module works out:

  * a per-project migration ``plan`` — for each project, whether it will be
    adopted in place or copied out, and to where; and
  * a list of structural ``issues`` (nested projects, destination collisions,
    case/unicode name clashes, shared/symlinked trees, unreadable sources, ...)
    that the UI surfaces and, where blocking, must have acknowledged.

It is deliberately free of Django and sqlite so it can be unit-tested against a
plain list of project dicts and a temp directory tree. Both ``SQLiteValidator``
and ``SQLiteImporter`` call :func:`analyse_structure` so validate and import can
never disagree about what will happen.

See docs/LEGACY_MIGRATION_WIRE_CONTRACT.md for the ``Issue`` / ``PlanEntry`` shapes.
"""

import os
import shutil
import unicodedata
from pathlib import Path

# Issue type constants (must match the wire contract enum).
NESTED_PROJECT = "nested_project"
DEST_COLLISION = "dest_collision"
SHARED_DIRECTORY = "shared_directory"
NOT_PROJECT_ROOT = "not_project_root"
DEST_UNWRITABLE = "dest_unwritable"
INSUFFICIENT_SPACE = "insufficient_space"
I1_DIR_OVERLAP = "i1_dir_overlap"
CASE_CLASH = "case_clash"
UNICODE_CLASH = "unicode_clash"
SHARED_SUBTREE = "shared_subtree"
UNREADABLE_SOURCE = "unreadable_source"

SEVERITY_BLOCKING = "blocking"
SEVERITY_WARNING = "warning"
SEVERITY_INFO = "info"

MODE_IN_PLACE = "in_place"
MODE_COPY = "copy"


def _norm_name(name):
    """Case- and unicode-fold a name for collision detection.

    macOS/Windows filesystems are case-insensitive, and macOS stores names in
    NFD, so two legacy names that differ only by case or normalisation form
    collide at the destination. Fold both dimensions for comparison only.
    """
    return unicodedata.normalize("NFC", (name or "")).casefold()


def _resolve(path):
    """Canonicalise a path (expanduser + resolve symlinks/.. ) without touching disk hard.

    ``Path.resolve(strict=False)`` collapses ``..`` and follows existing symlinks
    but tolerates missing tails, which is what we want for pre-flight.
    """
    if not path:
        return None
    try:
        return Path(path).expanduser().resolve(strict=False)
    except (OSError, RuntimeError, ValueError):
        return None


def _is_under(inner, outer):
    """True if ``inner`` is at or below ``outer`` (both already resolved)."""
    if inner is None or outer is None:
        return False
    if inner == outer:
        return False  # same dir is a duplicate, not a nesting — handled elsewhere
    try:
        return inner.is_relative_to(outer)
    except AttributeError:  # Python < 3.9 safety, though repo targets 3.9+
        try:
            inner.relative_to(outer)
            return True
        except ValueError:
            return False


def _unique_dest(root, name, taken):
    """Pick a collision-free destination dir under ``root`` for ``name``.

    Avoids clashes both with names already chosen this run (``taken``, a set of
    case/NFC-folded names) and with directories already present on disk. Mirrors
    the ``_legacy`` / ``_N`` suffixing used for project *names* elsewhere.
    """
    base = name or "project"
    candidates = [base, f"{base}_legacy"] + [f"{base}_legacy_{i}" for i in range(2, 1000)]
    for cand in candidates:
        folded = _norm_name(cand)
        if folded in taken:
            continue
        dest = root / cand
        if dest.exists():
            continue
        taken.add(folded)
        return dest, (cand if cand != base else None)
    # Pathological fallback — should never happen.
    taken.add(_norm_name(base))
    return root / base, None


def _looks_like_project_root(directory):
    """Cheap check that a directory is a CCP4i2 project root.

    A real legacy project has a ``CCP4_JOBS`` subdirectory. Missing it (when the
    directory otherwise exists) usually means a wrong/moved remap target.
    """
    if not directory:
        return True  # empty dir is a separate (dir_empty) concern, not ours
    p = Path(directory)
    if not p.is_dir():
        return True  # existence is checked separately; don't double-report
    return (p / "CCP4_JOBS").is_dir()


def _dir_readable(directory):
    """True if ``directory`` and a shallow sample of its tree are readable.

    A full walk of a large project is wasteful; sample the root plus the two
    content subdirs, which is where copytree would fail first.
    """
    p = Path(directory)
    for candidate in (p, p / "CCP4_JOBS", p / "CCP4_IMPORTED_FILES"):
        if candidate.exists() and not os.access(candidate, os.R_OK):
            return False
    return True


def _tree_size(directory, exclude=()):
    """Sum file sizes under ``directory``, skipping any path under ``exclude``."""
    total = 0
    exclude = [e for e in exclude if e is not None]
    root = Path(directory)
    for dirpath, _dirnames, filenames in os.walk(root, followlinks=False):
        here = _resolve(dirpath)
        if here is not None and any(_is_under(here, ex) or here == ex for ex in exclude):
            continue
        for fn in filenames:
            fp = Path(dirpath) / fn
            try:
                if fp.is_symlink():
                    continue
                total += fp.stat().st_size
            except OSError:
                pass
    return total


def analyse_structure(projects, *, copy_files, dest_root, check_disk=True):
    """Analyse legacy projects for structural migration hazards.

    Args:
        projects: list of dicts, each ``{"id", "name", "directory", "i1_directory"}``
            where ``directory`` / ``i1_directory`` are already remapped.
        copy_files: whether the run copies (True) or adopts in place (False).
        dest_root: ``Path`` (or str) — destination root for copied projects.
        check_disk: if False, skip filesystem probes (readability, space,
            project-root, tree walks) and reason purely about paths. Used by the
            plan-only path and by tests without a real tree.

    Returns:
        ``(issues, plan)`` where ``issues`` is a list of Issue dicts and ``plan``
        is a list of PlanEntry dicts (one per input project), both matching the
        wire contract.
    """
    dest_root = Path(dest_root).expanduser() if dest_root else None

    # Resolve directories once.
    resolved = {}       # project id -> resolved Path (or None)
    for p in projects:
        resolved[p["id"]] = _resolve(p.get("directory"))

    issues = []

    # ---- 1. Nesting: ordered pairs, inner under outer ------------------------
    # A project is "nested" (and must be copied out) if it lies under another.
    nested_ids = set()
    for a in projects:
        ra = resolved[a["id"]]
        if ra is None:
            continue
        for b in projects:
            if a["id"] == b["id"]:
                continue
            rb = resolved[b["id"]]
            if _is_under(rb, ra):
                # b is inside a  ->  copy the inner one (b) out
                nested_ids.add(b["id"])
                issues.append({
                    "type": NESTED_PROJECT,
                    "severity": SEVERITY_BLOCKING,
                    "mode_scope": "both",
                    "projects": [a["id"], b["id"]],
                    "detail": (
                        f"Project '{b['name']}' lies inside project '{a['name']}' "
                        f"({rb} under {ra})"
                    ),
                    # resolution filled in after the plan is built (needs dest)
                    "resolution": None,
                    "_inner": b["id"],
                })

    # ---- 3. Shared / duplicate resolved directories (past `unique`) ----------
    by_resolved = {}
    for pid, rp in resolved.items():
        if rp is not None:
            by_resolved.setdefault(rp, []).append(pid)
    for rp, pids in by_resolved.items():
        if len(pids) > 1:
            names = [pp["name"] for pp in projects if pp["id"] in pids]
            issues.append({
                "type": SHARED_DIRECTORY,
                "severity": SEVERITY_WARNING,
                "mode_scope": "both",
                "projects": pids,
                "detail": f"{len(pids)} projects resolve to the same directory {rp}: {', '.join(names)}",
                "resolution": None,
            })

    # ---- Build the per-project plan ------------------------------------------
    # A project is copied if the run is copy mode OR it is nested.
    plan = []
    taken_dest_names = set()   # case/NFC-folded dest basenames chosen this run
    dest_for = {}              # project id -> chosen dest Path (copy only)

    for p in projects:
        pid = p["id"]
        rp = resolved[pid]
        # Keep the ORIGINAL (remapped) directory string for the plan so an
        # in-place project's Project.directory is unchanged from legacy. The
        # resolved path (rp) is used only for nesting/collision reasoning.
        source_dir = p.get("directory") or ""
        source_exists = bool(rp and rp.is_dir()) if check_disk else bool(source_dir)
        is_nested = pid in nested_ids
        will_copy = bool(copy_files or is_nested)

        if will_copy and dest_root is not None:
            dest_path, renamed = _unique_dest(dest_root, p["name"], taken_dest_names)
            dest_for[pid] = dest_path
            reason = "nested" if (is_nested and not copy_files) else "copy_files"
            plan.append({
                "legacy_project_id": pid,
                "name": p["name"],
                "source_dir": source_dir,
                "mode": MODE_COPY,
                "dest": str(dest_path),
                "reason": reason,
                "source_exists": source_exists,
                "renamed_to": renamed,
            })
        else:
            plan.append({
                "legacy_project_id": pid,
                "name": p["name"],
                "source_dir": source_dir,
                "mode": MODE_IN_PLACE,
                "dest": source_dir,
                "reason": MODE_IN_PLACE,
                "source_exists": source_exists,
                "renamed_to": None,
            })

    # Fill nesting resolutions now that dests are known.
    for issue in issues:
        if issue["type"] == NESTED_PROJECT:
            inner = issue.pop("_inner")
            dest = dest_for.get(inner)
            if dest is not None:
                inner_name = next((pp["name"] for pp in projects if pp["id"] == inner), inner)
                issue["resolution"] = f"'{inner_name}' will be copied to {dest} and repointed"

    # ---- 2 / 7 / 8. Destination collisions (case + unicode aware) ------------
    # _unique_dest already prevents collisions among THIS run's copies, but a
    # chosen dest may still clash with a pre-existing directory, and two source
    # projects may differ only by case/normalisation. Report both.
    if copy_files or nested_ids:
        seen_folded = {}   # folded original name -> first project name
        for p in projects:
            entry = next(e for e in plan if e["legacy_project_id"] == p["id"])
            if entry["mode"] != MODE_COPY:
                continue
            folded = _norm_name(p["name"])
            if folded in seen_folded and seen_folded[folded] != p["name"]:
                clash_type = CASE_CLASH if p["name"].casefold() == seen_folded[folded].casefold() \
                    and p["name"] != seen_folded[folded] else UNICODE_CLASH
                issues.append({
                    "type": clash_type,
                    "severity": SEVERITY_BLOCKING,
                    "mode_scope": "copy",
                    "projects": [p["id"]],
                    "detail": (
                        f"Destination name '{p['name']}' clashes with '{seen_folded[folded]}' "
                        f"on a case/normalisation-insensitive filesystem "
                        f"(resolved by copying to {entry['dest']})"
                    ),
                    "resolution": f"copied to {entry['dest']}",
                })
            else:
                seen_folded.setdefault(folded, p["name"])

    if check_disk:
        # ---- 4. Not a project root -------------------------------------------
        for p in projects:
            d = p.get("directory")
            if d and Path(d).is_dir() and not _looks_like_project_root(d):
                issues.append({
                    "type": NOT_PROJECT_ROOT,
                    "severity": SEVERITY_WARNING,
                    "mode_scope": "both",
                    "projects": [p["id"]],
                    "detail": f"'{p['name']}' directory {d} has no CCP4_JOBS subdirectory",
                    "resolution": None,
                })

        # ---- 10. Unreadable source (copy only) -------------------------------
        for entry in plan:
            if entry["mode"] == MODE_COPY and entry["source_exists"]:
                if not _dir_readable(entry["source_dir"]):
                    issues.append({
                        "type": UNREADABLE_SOURCE,
                        "severity": SEVERITY_BLOCKING,
                        "mode_scope": "copy",
                        "projects": [entry["legacy_project_id"]],
                        "detail": f"'{entry['name']}' source {entry['source_dir']} is not readable",
                        "resolution": None,
                    })

        # ---- 9. Shared subtree via symlink (copy only) -----------------------
        # Two different projects whose CCP4_JOBS/CCP4_IMPORTED_FILES resolve into
        # a common real directory. Blocking only if the shared target lies under
        # another migrating project's dest.
        subtree_targets = {}   # resolved subtree -> [project ids]
        for entry in plan:
            if entry["mode"] != MODE_COPY or not entry["source_exists"]:
                continue
            for sub in ("CCP4_JOBS", "CCP4_IMPORTED_FILES"):
                rt = _resolve(Path(entry["source_dir"]) / sub)
                if rt is not None and rt.is_dir():
                    subtree_targets.setdefault(rt, []).append(entry["legacy_project_id"])
        for rt, pids in subtree_targets.items():
            if len(set(pids)) > 1:
                issues.append({
                    "type": SHARED_SUBTREE,
                    "severity": SEVERITY_WARNING,
                    "mode_scope": "copy",
                    "projects": list(set(pids)),
                    "detail": f"{len(set(pids))} projects share a content directory {rt} (symlinked)",
                    "resolution": None,
                })

        # ---- 5. Destination writable + free space (copy only) ----------------
        if (copy_files or nested_ids) and dest_root is not None:
            try:
                dest_root.mkdir(parents=True, exist_ok=True)
            except OSError:
                pass
            if not (dest_root.is_dir() and os.access(dest_root, os.W_OK)):
                issues.append({
                    "type": DEST_UNWRITABLE,
                    "severity": SEVERITY_BLOCKING,
                    "mode_scope": "copy",
                    "projects": [],
                    "detail": f"Destination root {dest_root} is not writable",
                    "resolution": None,
                })
            else:
                # Estimate total to copy (excluding nested inner trees from outer sums).
                nested_resolved = [resolved[i] for i in nested_ids if resolved.get(i)]
                needed = 0
                for entry in plan:
                    if entry["mode"] == MODE_COPY and entry["source_exists"]:
                        needed += _tree_size(entry["source_dir"], exclude=nested_resolved)
                try:
                    free = shutil.disk_usage(dest_root).free
                    if needed > free:
                        issues.append({
                            "type": INSUFFICIENT_SPACE,
                            "severity": SEVERITY_BLOCKING,
                            "mode_scope": "copy",
                            "projects": [],
                            "detail": (
                                f"Copy needs ~{needed // (1024*1024)} MiB but only "
                                f"{free // (1024*1024)} MiB free at {dest_root}"
                            ),
                            "resolution": None,
                        })
                except OSError:
                    pass

    # ---- 6. i1 directory overlap (warn only) ---------------------------------
    for p in projects:
        ri1 = _resolve(p.get("i1_directory"))
        if ri1 is None:
            continue
        for other in projects:
            if other["id"] == p["id"]:
                continue
            ro = resolved[other["id"]]
            if _is_under(ri1, ro) or ri1 == ro:
                issues.append({
                    "type": I1_DIR_OVERLAP,
                    "severity": SEVERITY_WARNING,
                    "mode_scope": "both",
                    "projects": [p["id"], other["id"]],
                    "detail": (
                        f"'{p['name']}' i1 directory {ri1} overlaps project "
                        f"'{other['name']}'"
                    ),
                    "resolution": None,
                })

    return issues, plan


def blocking_unresolved(issues):
    """Count issues that must be acknowledged (blocking, no auto-resolution)."""
    return sum(
        1 for i in issues
        if i["severity"] == SEVERITY_BLOCKING and not i.get("resolution")
    )


def nested_excludes_for(plan, resolved_by_id):
    """Map each copied project's id -> set of resolved inner dirs to exclude.

    When copying an outer project's tree, any *other* migrating project directory
    that falls under it must be excluded so its content is not duplicated (the
    inner project is copied separately to its own dest).
    """
    excludes = {}
    for entry in plan:
        if entry["mode"] != MODE_COPY:
            continue
        outer = _resolve(entry["source_dir"])
        if outer is None:
            continue
        inner_dirs = set()
        for other_id, other_resolved in resolved_by_id.items():
            if other_id == entry["legacy_project_id"] or other_resolved is None:
                continue
            if _is_under(other_resolved, outer):
                inner_dirs.add(other_resolved)
        excludes[entry["legacy_project_id"]] = inner_dirs
    return excludes
