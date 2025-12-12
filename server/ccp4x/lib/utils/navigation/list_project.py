import uuid
import os
from ccp4x.db import models


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
