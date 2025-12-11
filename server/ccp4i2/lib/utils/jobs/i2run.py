from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.fundamental_types import CList
from ccp4i2.core.base_object.cdata import CData
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4Data
from ccp4i2.core import CCP4Container

from ccp4x.db import models
from ..containers.get_container import get_job_container
from ..containers.find_objects import find_objects


def i2run_for_job(job: models.Job):
    container = get_job_container(job)
    if not container:
        return None
    command: str = f"{job.task_name} --project_name {job.project.name} "
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
            print(full_path, candidate_path)
            return candidate_path

    # If no unique shorter path found, return the full relative path
    print(full_path, ".".join(full_parts))
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
            print("grandchild", grandchild, element_text)
            if len(element_text) > 0:
                command += f" --{minimal_path(child.objectPath(), container)} "
                command += handle_element(grandchild)
        return command

    def handle_nonleaf_child(command, child, container):
        element_text = handle_element(child)
        if len(element_text) > 0:
            command += f" --{minimal_path(child.objectPath(), container)} "
            command += handle_element(child)
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
        command += f' --{minimal_path(child.objectPath(), container)} "{str(child)}" '
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
                results.append(f'"{path}={value}" ')

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
    return (
        candidate_path == full_path
        or len(
            find_objects(
                container,
                lambda x: x.objectPath().endswith(f".{candidate_path}"),
                multiple=True,
            )
        )
        == 1
    )
