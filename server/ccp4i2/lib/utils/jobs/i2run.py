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
from ..parameters.argument_names import argument_names

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


def _table_key(full_path: str, container: CCP4Container) -> str:
    """*full_path* rooted the way the argument table is rooted.

    The table is keyed from the container's own name down
    ("container.inputData.XYZIN"), but ``objectPath()`` is rooted at whatever
    is above it -- "freerflag.container.inputData.XYZIN" once the plugin that
    owns the container is in the picture, which it now always is, because
    get_job_container keeps the plugin alive on purpose.

    The old minimiser never noticed: it compared suffixes, so a differing root
    was invisible to it. An exact lookup has to be told.
    """
    prefix = container.objectPath()
    root = container.objectName()
    if full_path == prefix:
        return root
    if full_path.startswith(f"{prefix}."):
        return f"{root}.{full_path[len(prefix) + 1:]}"
    return full_path


def minimal_path(full_path, container: CCP4Container, names=None) -> str:
    """The name i2run accepts for the parameter at *full_path*.

    Looks the answer up in the one table that decides it
    (``lib.utils.parameters.argument_names``), which is also what i2run's
    argparse arguments are built from --- so a rendered command cannot name a
    parameter in a spelling the parser rejects.

    This used to re-derive the answer by walking the container and testing
    ``objectPath().endswith("." + candidate)`` for every candidate suffix of
    every parameter: a second implementation of the same rule, and a few
    hundred thousand predicate calls per render on a task the size of
    ``servalcat_pipe``.

    *names* is the precomputed table; it is built once per render and threaded
    through, because building it instantiates nothing but does walk the
    container.
    """
    if names is None:
        names = argument_names(container)

    name = names.get(full_path) or names.get(_table_key(full_path, container))
    if name is not None:
        return name

    # Not a parameter i2run exposes (it should be, for anything we render).
    # Fall back to the path relative to the container rather than inventing a
    # spelling, and say so, because this is the shape of a real bug.
    logger.warning(
        "No i2run argument name for %s; falling back to the relative path",
        full_path,
    )
    parts = full_path.split(".")
    return ".".join(parts[1:]) if len(parts) > 1 else full_path


def _is_list(object: CData) -> bool:
    return isinstance(object, (CList, CCP4Data.CList))


def _is_container(object: CData) -> bool:
    return isinstance(object, (CCP4Container.CContainer, CContainer))


def _is_file(object: CData) -> bool:
    return isinstance(object, (CCP4File.CDataFile, CDataFile))


def _file_is_registered(node) -> bool:
    """Is this file identified by a database id?"""
    try:
        db_id = getattr(node, "dbFileId", None)
        return db_id is not None and db_id.isSet()
    except Exception:
        return False


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
    command: str,
    element: CData,
    container: CCP4Container,
    exclude: list[str] = None,
    names: dict = None,
) -> str:
    if exclude is None:
        exclude = []
    if names is None:
        names = argument_names(container)

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
                command += f" --{minimal_path(child.objectPath(), container, names)}"
                command += f" {element_text}"
        return command

    def handle_nonleaf_child(command, child, container):
        element_text = handle_element(child)
        if len(element_text) > 0:
            command += f" --{minimal_path(child.objectPath(), container, names)}"
            command += f" {element_text}"
        return command

    def process_child(command, child, container, exclude):
        if should_skip_child(child, exclude):
            return command
        if is_unset_nonlist_noncontainer(child):
            return command
        if _is_container(child):
            return extend_i2run(command, child, container, exclude, names)
        if _is_list(child):
            return handle_list_child(command, child, container)
        if not _is_leaf(child):
            return handle_nonleaf_child(command, child, container)
        command += f' --{minimal_path(child.objectPath(), container, names)} "{str(child)}"'
        return command

    for child in element.children():
        command = process_child(command, child, container, exclude)

    return command


def _file_use_text(node) -> str:
    """``fileUse=[N].PARAM`` for a file another job produced, else "".

    This is what makes a rendered command editable. The first thing anyone does
    with a surfaced i2run call is change its inputs, and nobody can retype
    ``dbFileId=a3ed78ad466845a88765271c38d15149`` or work out what it was --
    whereas ``[3].XYZOUT`` says which job and which output, and edits cleanly
    to ``[-1].XYZOUT`` or ``prosmart_refmac[-1].XYZOUT`` for a script.

    An ABSOLUTE reference is rendered on purpose. Relative ones are for people
    to write: a command that meant "the latest refmac" would quietly resolve
    somewhere else next week, and a surfaced command should reproduce the job
    it was surfaced from.

    Empty for an imported file, which has no producing job to name.
    """
    from ....db import models
    from ..files.file_use import file_use_for_file

    try:
        db_id = str(node.dbFileId)
    except Exception:
        return ""
    if not db_id:
        return ""

    the_file = models.File.objects.filter(uuid=db_id).select_related("job").first()
    if the_file is None:
        return ""
    reference = file_use_for_file(the_file)
    return f'"fileUse={reference}"' if reference else ""


def handle_element(item: CData) -> str:
    # If this is a simple element, then simply return the corresponding quoted string value
    if _is_leaf(item):
        return f'"{str(item)}"'

    if _is_file(item) and _file_is_registered(item):
        file_use = _file_use_text(item)
        if file_use:
            return file_use

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

            # A registered file needs its dbFileId and nothing else. The id
            # identifies the row, and the row carries the rest -- so emitting
            # project/baseName/relPath alongside it is redundant, and worse
            # than redundant: it makes an inconsistent command representable
            # (edit baseName, leave dbFileId, and the two now disagree about
            # which file is meant).
            #
            # contentFlag is deliberately omitted too, and that one matters.
            # It is set by INTROSPECTION when the file is registered, so a
            # command restating it can contradict the file it names -- asking
            # anyone to supply it by hand for a database file is dangerous.
            # Verified: `--F_SIGF "dbFileId=<uuid>"` alone runs a job to
            # completion, contentFlag and all.
            if _is_file(node) and _file_is_registered(node):
                child_nodes = [
                    child
                    for child in child_nodes
                    if child.objectName() == "dbFileId"
                ]

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
