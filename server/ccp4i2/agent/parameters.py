"""A job's parameters as an agent can use them: one line each.

The container's JSON (the ``container`` endpoint) is complete and large —
some 300 KB for a Phaser job, most of it the qualifiers of every sub-field
of every file. An agent needs, per parameter: where it is, what the interface
calls it, its value, whether it must be set, and its choices.
"""
from ..core.base_object.cdata_file import CDataFile
from ..core.base_object.fundamental_types import CList
from ..core.base_object.ccontainer import CContainer

SECTIONS = ("inputData", "controlParameters", "outputData")
_BASIC = (str, int, float, bool, type(None))


def _qualifier(obj, key):
    try:
        value = obj.get_qualifier(key)
    except Exception:  # a class with an unusual qualifier store
        return None
    return value if value not in ("", [], NotImplemented) else None


def _plain(value):
    if isinstance(value, _BASIC):
        return value
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return str(value)


def _is_set(obj):
    try:
        return bool(obj.isSet())
    except Exception:
        return False


def _value(obj):
    if isinstance(obj, CDataFile):
        if not _is_set(obj):
            return None
        annotation = getattr(obj, "annotation", None)
        out = {"file": str(getattr(obj, "baseName", "") or "")}
        if annotation is not None and str(annotation):
            out["annotation"] = str(annotation)
        db_file_id = getattr(obj, "dbFileId", None)
        if db_file_id is not None and str(db_file_id):
            out["fileId"] = str(db_file_id)
        return out
    if isinstance(obj, CList):
        return [_value(item) for item in obj]
    children = [c for c in obj.children() if c.objectName() and not c.objectName().startswith("_")]
    if children:  # a composite (an ensemble, a cell): its set fields
        fields = {c.objectName(): _value(c) for c in children if _is_set(c)}
        return {k: v for k, v in fields.items() if v not in (None, [], "")} or None
    try:
        value = obj.get()
    except Exception:
        value = str(obj)
    return _plain(value)


def describe(obj, path):
    entry = {"path": path, "class": type(obj).__name__}
    label = _qualifier(obj, "guiLabel")
    if label:
        entry["label"] = str(label)
    tip = _qualifier(obj, "toolTip")
    if tip:
        entry["tip"] = str(tip)
    entry["set"] = _is_set(obj)
    value = _value(obj)
    if value not in (None, [], ""):
        entry["value"] = value
    if _qualifier(obj, "allowUndefined") is False:
        entry["required"] = True
    enumerators = _qualifier(obj, "enumerators")
    if enumerators:
        entry["choices"] = _plain(enumerators)
        menu = _qualifier(obj, "menuText")
        if menu and len(menu) == len(enumerators):
            entry["choice_labels"] = _plain(menu)
    default = _qualifier(obj, "default")
    if default is not None and not isinstance(obj, CDataFile):
        entry["default"] = _plain(default)
    if isinstance(obj, CDataFile):
        mime = _qualifier(obj, "mimeTypeName")
        if mime:
            entry["file_type"] = str(mime)
    return entry


def _walk(container, prefix):
    for child in container.children():
        name = child.objectName()
        if not name or name.startswith("_"):
            continue
        path = f"{prefix}.{name}"
        if isinstance(child, CContainer):
            yield from _walk(child, path)
        else:
            yield child, path


def summarise(container, sections=SECTIONS, only_set=False, query=None):
    """The parameters of a job's container, one dict each, in def.xml order.

    ``sections`` limits it to some of inputData / controlParameters /
    outputData; ``only_set`` to those with a value; ``query`` to those whose
    path or label contains it (case-insensitive).
    """
    query = query.lower() if query else None
    out = []
    for section in container.children():
        name = section.objectName()
        if name not in sections or not isinstance(section, CContainer):
            continue
        for obj, path in _walk(section, name):
            try:
                entry = describe(obj, path)
            except Exception as err:  # never let one odd class hide the rest
                entry = {"path": path, "class": type(obj).__name__, "error": str(err)}
            if only_set and not entry.get("set"):
                continue
            if query and query not in path.lower() and query not in entry.get("label", "").lower():
                continue
            out.append(entry)
    return out
