"""Add one item to a list parameter of a pending job.

A list is otherwise set whole, as a JSON array of items shaped as its
item_fields say. For a nested list (a Phaser ensemble holds a list of
models) that is a structure a person builds in the interface by clicking
"add" and filling in each row; an agent asked to write it out whole failed
some forty times on BetaBlip (gpt-oss:20b, 2026-10-12) while trying, over
and over, to do what the interface does: add an item, then fill it in.

The item is the list's own (``CList.addItem``, so its class, qualifiers and
defaults are the list's), optionally with some of its fields set, and the
answer names the paths of what it holds, so the next step needs no guessing.
"""
import logging
from typing import Any, Dict, Optional

from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core.base_object.fundamental_types import CList
from ccp4i2.db import models
from ccp4i2.lib.response import Result
from ccp4i2.lib.utils.plugins.plugin_context import get_plugin_with_context

from .set_param import normalize_object_path

logger = logging.getLogger(__name__)


def item_paths(item, path: str) -> Dict[str, list]:
    """What an item holds, by path: its plain fields, its file fields and
    the lists inside it. An item that is itself a file is its own file."""
    if isinstance(item, CDataFile):
        return {"fields": [], "files": [path], "lists": []}
    out = {"fields": [], "files": [], "lists": []}
    for child in item.children():
        name = child.objectName()
        if not name or name.startswith("_"):
            continue
        here = f"{path}.{name}"
        if isinstance(child, CList):
            out["lists"].append(here)
        elif isinstance(child, CDataFile):
            out["files"].append(here)
        else:
            out["fields"].append(here)
    return out


def add_list_item(job: models.Job, object_path: str,
                  values: Optional[Dict[str, Any]] = None) -> Result[Dict[str, Any]]:
    """Append an item to the list at ``object_path`` of a pending job, set
    any of its fields given in ``values`` ({field: value}, a file as the
    set_parameter endpoint takes one), and save. Nothing is kept if a value
    is refused."""
    if job.status not in [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]:
        return Result.fail(
            f"Cannot modify parameters on job with status '{job.status}'. "
            f"Use clone API to create an editable copy first.")
    plugin_result = get_plugin_with_context(job)
    if not plugin_result.success:
        return Result.fail(f"Failed to load plugin: {plugin_result.error}",
                           details=plugin_result.error_details)
    plugin = plugin_result.data
    normalized = normalize_object_path(object_path, job.task_name)
    try:
        target = plugin.container.find_by_path(normalized, skip_first=True)
    except AttributeError:
        target = None
    if not isinstance(target, CList):
        what = "nothing" if target is None else f"a {type(target).__name__}, not a list"
        return Result.fail(f"'{object_path}' is {what}: add_list_item adds to a list parameter")

    item = target.addItem()
    index = len(target) - 1
    # The path in the caller's section form (inputData.ENSEMBLES[1]), whatever form came in
    section_path = normalized.split(".", 1)[1] if normalized.startswith(f"{job.task_name}.") else normalized
    item_path = f"{section_path}[{index}]"
    if values:
        from .set_parameter import file_type_errors, shape_errors
        for field, value in values.items():
            field_path = f"{normalized}[{index}].{field}"
            try:
                child = plugin.container.find_by_path(field_path, skip_first=True)
                problems = (shape_errors(child, value, f"{item_path}.{field}")
                            + file_type_errors(child, value, f"{item_path}.{field}"))
                if problems:
                    raise ValueError("; ".join(problems))
                set_obj = plugin.container.set_parameter(field_path, value, skip_first=True)
                from .set_parameter import after_api_set
                after_api_set(set_obj)
            except (AttributeError, ValueError, TypeError) as err:
                target.pop()
                fields = item_paths(item, item_path)
                names = [p.rsplit(".", 1)[-1] for p in fields["fields"] + fields["files"] + fields["lists"]]
                return Result.fail(f"Item not added: field '{field}' of {item_path}: {err}. "
                                   f"Its fields: {', '.join(names)}")

    error = plugin.saveDataToXml(str(job.directory / "input_params.xml"))
    if error and hasattr(error, "hasError") and error.hasError():
        logger.error("Failed to save parameters for job %s: %s", job.uuid, error)
    logger.info("Added item %s to job %s", item_path, job.uuid)
    return Result.ok({"item": item_path, "index": index, "length": len(target),
                      **item_paths(item, item_path)})
