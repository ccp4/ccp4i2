import logging
import os
import shlex
import sys
from pathlib import Path

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.fundamental_types import CList
from ccp4i2.core.base_object.cdata import CData
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4Data
from ccp4i2.core import CCP4Container

from ccp4i2.db import models
from ..containers.get_container import get_job_container

logger = logging.getLogger(f"ccp4i2:{__name__}")


def i2run_for_job(job: models.Job):
    container = get_job_container(job)
    # `is None`, not falsiness: a CData's truth value is its own business (an
    # empty container and a CBoolean(False) are both falsy for reasons that
    # have nothing to do with whether we were handed one).
    if container is None:
        return None
    command: str = f"{job.task_name} --project_name {shlex.quote(job.project.name)}"
    command = extend_i2run(
        command,
        container,
        container,
        exclude=[
            "outputData",
            "guiAdmin",
            "guiControls",
            "patchSelection",
            "guiParameters",
            "temporary",
        ],
    )
    return command


def ccp4i2_root() -> Path:
    """The installed ``ccp4i2`` package directory."""
    import ccp4i2

    return Path(ccp4i2.__file__).parent


def i2run_working_directory():
    """The directory the rendered command must be run from, or None.

    In a pip-installed tree, ``-m ccp4i2.cli.i2run`` resolves wherever you
    are. In a dev checkout it does not: CCP4's own bundle ships a legacy
    ``ccp4i2`` directory with no ``__init__.py``, i.e. a namespace-package
    portion of the same name, and from an unrelated working directory an
    editable install loses the race for ``ccp4i2.core``. Running from the
    directory holding ``manage.py`` puts the checkout first on ``sys.path``
    and settles it.
    """
    server_dir = ccp4i2_root().parent
    if (server_dir / "manage.py").is_file():
        return server_dir
    return None


# Variables the desktop app sets for the server process that change where
# i2run looks or how it behaves, and that a fresh terminal will not have.
# Deliberately not the server's own plumbing (UVICORN_PORT, NEXT_ADDRESS, the
# parent-pid watchdog, MPLCONFIGDIR): those belong to the running server, not
# to a job.
I2RUN_ENVIRONMENT_VARIABLES = (
    # Where the database and projects live. If the GUI is running with either
    # of these and the terminal is not, the command silently addresses a
    # DIFFERENT database -- the one failure mode worth going out of our way to
    # prevent, because it looks like it worked.
    "CCP4I2_HOME",
    "CCP4I2_PROJECTS_DIR",
    "CCP4I2_DB_FILE",
    # Where jobs run. The desktop always sets this to local explicitly.
    "CCP4I2_JOB_TARGET",
    # Keeps matplotlib (report graphs) off a GUI backend.
    "MPLBACKEND",
)


def i2run_environment() -> dict:
    """The variables this server runs with that a terminal would need too.

    Only what is actually set: reporting a default as though it were a setting
    would be noise, and worse, would go stale the day the default changes.
    """
    return {
        name: os.environ[name]
        for name in I2RUN_ENVIRONMENT_VARIABLES
        if os.environ.get(name)
    }


def i2run_ccp4_setup():
    """Path to the CCP4 setup script to source first, or None.

    None on Windows, which has no ``ccp4.setup-sh``: there the desktop app
    builds the CCP4 environment itself (client/main/ccp4i2-setup-windows.ts),
    and a person's equivalent is the CCP4 command prompt.
    """
    ccp4 = os.environ.get("CCP4")
    if not ccp4:
        return None
    script = Path(ccp4) / "bin" / "ccp4.setup-sh"
    return str(script) if script.is_file() else None


def i2run_command_line(job: models.Job):
    """Render *job* as a command someone can paste into a terminal.

    Returns ``(working_directory, arguments, command_line)``: the directory to
    run from (None when it runs from anywhere), the argument list i2run itself
    parses, and the whole invocation.

    ``-m ccp4i2.cli.i2run`` rather than the ``i2run`` console script, which is
    shadowed on any normal CCP4 setup by the legacy Qt ``$CCP4/bin/i2run``,
    and rather than ``manage.py i2run``, which only exists in a checkout.
    """
    arguments = i2run_for_job(job)
    if not arguments:
        return None, None, None

    command_line = f"ccp4-python -m ccp4i2.cli.i2run {arguments}"
    return i2run_working_directory(), arguments, command_line


def minimal_path(full_path, container: CCP4Container) -> str:
    """
    Get the minimal unique path of a container relative to another container.
    Starts with the last element and adds path elements until the path is unique.
    """
    full_parts = full_path.split(".")

    # Start with the last element and gradually add more elements
    for i in range(1, len(full_parts) + 1):
        # Take the last i elements
        candidate_path_parts = full_parts[-i:]
        candidate_path = ".".join(candidate_path_parts)

        # Test if this path is unique within the container
        if _is_path_unique(candidate_path, container, full_path):
            logger.debug("%s -> %s", full_path, candidate_path)
            return candidate_path

    # If no unique shorter path found, return the full relative path
    logger.debug("%s -> %s (not shortened)", full_path, ".".join(full_parts))
    return ".".join(full_parts)


def _is_list(object: CData) -> bool:
    return isinstance(object, (CList, CCP4Data.CList))


def _is_container(object: CData) -> bool:
    return isinstance(object, (CCP4Container.CContainer, CContainer))


def _is_file(object: CData) -> bool:
    return isinstance(object, (CCP4File.CDataFile, CDataFile))


def _is_leaf(object: CData) -> bool:
    """
    Return True if the object is a leaf node.

    A leaf node is one that either:
    1. Has no children() method, or
    2. has a children() method that returns length zero
    """
    # Check if object has children method
    if not hasattr(object, "children"):
        return True

    try:
        children = object.children()

        # If no children, it's a leaf
        if not children or len(children) == 0:
            return True

    except (AttributeError, TypeError):
        # If children() method fails, treat as leaf
        return True

    return False


def extend_i2run(
    command: str, element: CData, container: CCP4Container, exclude: list[str] = None
) -> str:
    if exclude is None:
        exclude = []

    def should_skip_child(child, exclude):
        return child.objectName() in exclude or child.objectName() == "temporary"

    def is_unset_nonlist_noncontainer(child):
        return (
            hasattr(child, "isSet")
            and not _is_list(child)
            and not _is_container(child)
            and not child.isSet(allowDefault=False)
        )

    def handle_list_child(command, child, container):
        for grandchild in child:
            element_text = handle_element(grandchild)
            if len(element_text) > 0:
                command += f" --{minimal_path(child.objectPath(), container)}"
                command += f" {element_text}"
        return command

    def handle_nonleaf_child(command, child, container):
        element_text = handle_element(child)
        if len(element_text) > 0:
            command += f" --{minimal_path(child.objectPath(), container)}"
            command += f" {element_text}"
        return command

    def process_child(command, child, container, exclude):
        if should_skip_child(child, exclude):
            return command
        if is_unset_nonlist_noncontainer(child):
            return command
        if _is_container(child):
            return extend_i2run(command, child, container, exclude)
        if _is_list(child):
            return handle_list_child(command, child, container)
        if not _is_leaf(child):
            return handle_nonleaf_child(command, child, container)
        command += f' --{minimal_path(child.objectPath(), container)} "{str(child)}"'
        return command

    for child in element.children():
        command = process_child(command, child, container, exclude)

    return command


def handle_element(item: CData) -> str:
    # If this is a simple element, then simply return the corresponding quoted string value
    if _is_leaf(item):
        return f'"{str(item)}"'

    def traverse(node, path_parts, is_root=False):
        results = []
        # Don't include root node's objectName in path_parts
        next_path_parts = path_parts if is_root else path_parts + [node.objectName()]
        if _is_leaf(node):
            path = "/".join(next_path_parts)
            value = str(node)
            if (
                len(value) > 0
                and node is not None
                and hasattr(node, "_value")
                and node._value is not None
            ):
                results.append(f'"{path}={value}"')

        elif _is_list(node):
            for i_item, item in enumerate(node):
                next_path_parts[-1] = next_path_parts[-1] + f"[{i_item}]"
                results.extend(traverse(item, next_path_parts, is_root=True))

        else:
            child_nodes = node.children() if hasattr(node, "children") else []
            filtered_nodes = [
                child
                for child in child_nodes
                if not (
                    (
                        _is_file(node)
                        and child.objectName()
                        in ["fileContent", "subType", "annotation"]
                    )
                    or callable(child)  # Filter out executable/function children
                )
            ]
            for child in filtered_nodes:
                results.extend(traverse(child, next_path_parts, is_root=False))
        return results

    leaf_texts = traverse(item, [], is_root=True)
    return " ".join(leaf_texts)


def _is_path_unique(candidate_path, container, full_path):
    """Does *candidate_path* name exactly one object in *container*?

    Uses ``find_children_matching``, the traversal CData actually offers. The
    previous implementation went through ``lib.utils.containers.find_objects``,
    which reads ``within.CONTENTS`` --- an attribute modern ``CContainer`` does
    not have, so every call raised ``AttributeError`` and the whole
    i2run-command endpoint returned a 400.
    """
    if candidate_path == full_path:
        return True

    suffix = f".{candidate_path}"
    matches = container.find_children_matching(
        lambda x: hasattr(x, "objectPath") and x.objectPath().endswith(suffix)
    )
    return len(matches) == 1
