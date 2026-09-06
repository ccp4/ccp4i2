"""Removing a project's directory from disk - only on request, and only when it
is recognisably a project directory.

``DELETE /projects/<id>/`` removes a project's database records and leaves its
directory exactly as it is (the Qt-era default, so the project can be brought
back with Import -> a project folder). Only ``?delete_files=true`` removes the
directory too, through :func:`remove_project_directory`, which refuses anything
that does not look like a CCP4i2 project: ``Project.directory`` is a plain
text column, and an ``rmtree`` on a stale or mistaken path must not be able to
take a home directory with it.
"""

import logging
import shutil
from pathlib import Path
from typing import Optional, Tuple

from django.conf import settings

logger = logging.getLogger(__name__)

#: Any one of these marks a directory as a CCP4i2 project directory.
PROJECT_MARKERS = (
    "CCP4_JOBS",
    "CCP4_PROJECT_FILES",
    "CCP4_IMPORTED_FILES",
    "DATABASE.db.xml",
)


def project_directory_removal_problem(directory: Path) -> Optional[str]:
    """Why ``directory`` must not be removed, or ``None`` if it may be."""
    if not directory.is_absolute():
        return "path is not absolute"
    if directory.is_symlink():
        return "path is a symbolic link"
    if not directory.exists():
        return "directory does not exist"
    if not directory.is_dir():
        return "path is not a directory"
    resolved = directory.resolve()
    if resolved == Path(resolved.anchor):
        return "path is a filesystem root"
    if resolved == Path.home().resolve():
        return "path is the home directory"
    projects_dir = getattr(settings, "CCP4I2_PROJECTS_DIR", None)
    if projects_dir is not None and resolved == Path(projects_dir).resolve():
        return "path is the projects directory itself"
    if not any((directory / marker).exists() for marker in PROJECT_MARKERS):
        return (
            "directory does not look like a CCP4i2 project "
            f"(none of {', '.join(PROJECT_MARKERS)} present)"
        )
    return None


def remove_project_directory(directory: Path) -> Tuple[bool, Optional[str]]:
    """Remove ``directory`` if it passes the checks above.

    Returns ``(removed, reason)``: ``(True, None)`` when the directory is gone,
    ``(False, why)`` when it was left alone - including when ``rmtree`` itself
    failed part-way, in which case whatever remains is still on disk.
    """
    problem = project_directory_removal_problem(directory)
    if problem is not None:
        logger.warning("Not removing project directory %s: %s", directory, problem)
        return False, problem
    try:
        shutil.rmtree(directory)
    except OSError as exc:
        logger.warning("Failed to remove project directory %s: %s", directory, exc)
        return False, f"could not remove directory: {exc}"
    logger.warning("Removed project directory %s", directory)
    return True, None
