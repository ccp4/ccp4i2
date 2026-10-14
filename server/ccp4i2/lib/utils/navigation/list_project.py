import uuid
import os
from pathlib import Path

from ccp4i2.db import models


def get_directory_tree(
    path, max_depth=10, current_depth=0, max_files=10000, file_count=None
):
    if file_count is None:
        file_count = {"count": 0}

    if current_depth >= max_depth:
        return [{"error": f"Maximum depth ({max_depth}) exceeded"}]

    if file_count["count"] >= max_files:
        return [{"error": f"Maximum file count ({max_files}) exceeded"}]

    tree = []
    try:
        for entry in os.scandir(path):
            if file_count["count"] >= max_files:
                break

            try:
                stats = entry.stat(follow_symlinks=False)
                node = {
                    "path": entry.path,
                    "name": entry.name,
                    "type": "directory" if entry.is_dir() else "file",
                    "size": stats.st_size,
                    "mode": stats.st_mode,
                    "inode": stats.st_ino,
                    "device": stats.st_dev,
                    "nlink": stats.st_nlink,
                    "uid": stats.st_uid,
                    "gid": stats.st_gid,
                    "atime": stats.st_atime,
                    "mtime": stats.st_mtime,
                    "ctime": stats.st_ctime,
                }

                file_count["count"] += 1

                if entry.is_dir():
                    node["contents"] = get_directory_tree(
                        entry.path, max_depth, current_depth + 1, max_files, file_count
                    )
                tree.append(node)
            except (PermissionError, FileNotFoundError) as e:
                tree.append({"name": entry.name, "error": str(e)})
    except Exception as e:
        return [{"error": str(e)}]

    return tree


def contained_directory(root, relative=""):
    """The directory ``relative`` names under ``root``, refusing escapes.

    Resolved before the check, as ``project_file`` does, so neither ``..``
    nor a symlink pointing out of the tree can list what lies beyond it.
    """
    base = Path(root).resolve()
    target = (base / (relative or "")).resolve()
    if target != base and base not in target.parents:
        raise ValueError(f"{relative!r} is not inside the directory being browsed")
    return target


def list_directory_level(root, relative=""):
    """One directory's own entries, for a browser that opens folders lazily.

    The recursive walk stats every file of a tree on every poll: some
    thirteen thousand entries for a PanDDA run over 700 datasets, each a
    round trip when the projects live on a network share. A browser shows
    one folder at a time, so it asks for one folder at a time. Only files
    are stat'ed (for size and mtime); a directory's type comes from the
    directory entry itself, and whether it has children is left to the
    moment it is opened.

    A directory that does not exist yet (a job not started) lists as empty.
    """
    target = contained_directory(root, relative)
    if not target.is_dir():
        return []
    items = []
    for entry in os.scandir(target):
        try:
            is_dir = entry.is_dir()
        except OSError:
            is_dir = False
        node = {
            "path": entry.path,
            "name": entry.name,
            "type": "directory" if is_dir else "file",
        }
        if entry.is_symlink():
            node["link"] = True
        if not is_dir:
            try:
                stats = entry.stat(follow_symlinks=False)
                node["size"] = stats.st_size
                node["mtime"] = stats.st_mtime
            except OSError:
                pass
        items.append(node)
    return items


def search_directory(root, term, limit=500, max_entries=200000):
    """Entries under ``root`` whose names contain ``term``, case-insensitively.

    A browser that loads folders lazily cannot search what it has not
    loaded, so the search walks on the server, reading names only: no stat
    and no descent through symlinked directories. ``truncated`` says the
    answer stopped at ``limit`` matches or ``max_entries`` entries read.
    """
    base = contained_directory(root)
    needle = term.lower()
    matches = []
    seen = 0
    truncated = False
    pending = [base]
    while pending and not truncated:
        directory = pending.pop()
        try:
            entries = list(os.scandir(directory))
        except OSError:
            continue
        for entry in entries:
            seen += 1
            try:
                is_dir = entry.is_dir(follow_symlinks=False)
            except OSError:
                is_dir = False
            if needle in entry.name.lower():
                matches.append({
                    "path": entry.path,
                    "name": entry.name,
                    "type": "directory" if is_dir else "file",
                    "relative": os.path.relpath(entry.path, base),
                })
            if is_dir:
                pending.append(entry.path)
            if len(matches) >= limit or seen >= max_entries:
                truncated = True
                break
    matches.sort(key=lambda item: item["relative"])
    return {"items": matches, "truncated": truncated}


def browse(root, params):
    """Answer a lazy browser's request for ``root``, or None for the old walk.

    ``?path=<relative>`` lists that one folder; ``?search=<term>`` searches
    the tree by name. With neither, the caller serves the recursive walk,
    which stays for clients that still expect the whole tree.
    """
    if "search" in params:
        return search_directory(root, params["search"])
    if "path" in params:
        return list_directory_level(root, params["path"])
    return None


def job_directory(the_job_id):
    """A job's directory, by job id or uuid."""
    query = {"uuid": the_job_id} if _looks_like_uuid(the_job_id) else {"pk": the_job_id}
    return models.Job.objects.get(**query).directory


def list_job(the_job_id, max_depth=10, max_files=50000):
    """List one job's own directory, by job id or uuid.

    The project-wide listing walks every job in the project: on a campaign
    parent whose PanDDA job holds a directory per dataset that is thousands
    of entries to deliver a few, once per poll while the job runs. A job's
    tree is what the interface actually shows, so serve that.

    Returns the contents of the job directory itself (not a node wrapping
    it), so the caller reads it exactly as it reads the project listing's
    job node. A job whose directory does not exist yet lists as empty.
    """
    directory = job_directory(the_job_id)
    if not os.path.isdir(directory):
        return []
    return get_directory_tree(str(directory), max_depth=max_depth, max_files=max_files)


def _looks_like_uuid(value) -> bool:
    try:
        uuid.UUID(str(value))
    except (ValueError, AttributeError, TypeError):
        return False
    return True


def list_project(the_project_uuid: str, max_depth=10, max_files=50000):
    """
    List all files and directories in a project directory.

    Args:
        the_project_uuid: UUID of the project
        max_depth: Maximum depth to traverse (default: 10)
        max_files: Maximum number of files to process (default: 50000)

    Returns:
        Directory tree structure with files and subdirectories

    Note:
        If max_files limit is exceeded during traversal, some directories
        may not have their contents populated. Consider increasing the limit
        for large projects.
    """
    the_project = models.Project.objects.get(uuid=uuid.UUID(the_project_uuid))
    return get_directory_tree(
        str(the_project.directory), max_depth=max_depth, max_files=max_files
    )
