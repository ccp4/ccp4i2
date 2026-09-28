import uuid
import os
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
    query = {"uuid": the_job_id} if _looks_like_uuid(the_job_id) else {"pk": the_job_id}
    the_job = models.Job.objects.get(**query)
    directory = the_job.directory
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
