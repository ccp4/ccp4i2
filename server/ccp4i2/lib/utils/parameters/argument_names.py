"""How a task's parameters are named on an i2run command line.

A task's parameters live at dotted paths inside its container
(``container.inputData.XYZIN``), but nobody types those in full: i2run accepts
the shortest suffix that identifies a parameter uniquely, so ``--F_SIGF``
usually does, while ``servalcat_pipe`` genuinely needs
``--container.inputData.XYZIN`` because ``metalCoordWrapper.inputData.XYZIN``
exists too.

There is exactly one right answer to "what may I call this parameter", and this
module is it. It used to be computed twice: once here, by i2run's
``KeywordExtractor``, to build the argparse arguments, and once in
``lib/utils/jobs/i2run.py``, by walking the live container, to render a job as
a command. The two agreed --- 735 of 735 parameters across five tasks --- but
only because both happened to minimise against *every* leaf in the container,
including ``outputData`` and ``guiAdmin``. Nothing recorded that invariant, and
the day someone "tidied" one side to skip the output sections, it would have
shortened half the names on that side only, and the rendered command would have
died on an unrecognised argument.

The consumers are i2run's argparse builder (which needs the full keyword
dictionaries, with types and qualifiers) and the command renderer (which needs
only path -> name). Both get it from here.

This lives in ``lib`` rather than ``cli`` deliberately: ``cli/i2run`` already
imports ``lib.utils``, so cli -> lib is an existing edge, whereas lib -> cli
would be a new backwards one and a cycle waiting to happen.
"""

import logging
from typing import Any, Dict, List

from ccp4i2.core.CCP4Container import CContainer

logger = logging.getLogger(f"ccp4i2:{__name__}")


def leaf_paths(container: CContainer) -> List[Dict[str, Any]]:
    """Every parameter in *container*, as ``{path, object, qualifiers}``.

    Recursion stops at anything that is not a ``CContainer``, so a
    ``CDataFile`` or a ``CList`` is one parameter, not one per sub-field ---
    which is the granularity a command line names things at.

    Paths are rooted at the container's own name, so they read
    ``container.inputData.XYZIN`` and match ``CData.objectPath()``.
    """

    def traverse(node, path_parts):
        results = []
        if isinstance(node, CContainer):
            for child in node.children():
                results.extend(traverse(child, path_parts + [child.objectName()]))
        else:
            qualifiers = {}
            if hasattr(node, "get_merged_metadata"):
                meta = node.get_merged_metadata("qualifiers")
                if meta:
                    qualifiers = meta
            results.append(
                {"path": ".".join(path_parts), "object": node, "qualifiers": qualifiers}
            )
        return results

    all_leaves = traverse(container, [container.objectName()])

    # Deduplicate by path, keeping the first occurrence.
    unique: Dict[str, Dict[str, Any]] = {}
    for leaf in all_leaves:
        if leaf["path"] not in unique:
            unique[leaf["path"]] = leaf
    return list(unique.values())


def compute_minimum_paths(keywords: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Annotate *keywords* in place with the shortest name that identifies each.

    Adds:

    ``minimumPath``
        the shortest dot-separated suffix matching this parameter and no other,
        falling back to the full path when even that is not unique.
    ``simpleName``
        the last path element.
    ``isAmbiguousSimpleName`` / ``isShortestForSimpleName``
        whether the simple name is shared, and whether this is the shortest
        path bearing it --- i2run uses these to keep backward-compatible
        aliases for names that used to be unambiguous.

    Suffixes are compared element-wise, never as strings: a string ``endswith``
    can match halfway through a path element.
    """
    paths = [kw["path"].split(".") for kw in keywords]

    simple_name_map: Dict[str, List[int]] = {}
    for i, kw in enumerate(keywords):
        simple_name = paths[i][-1]
        kw["simpleName"] = simple_name
        simple_name_map.setdefault(simple_name, []).append(i)

    for indices in simple_name_map.values():
        if len(indices) > 1:
            shortest_idx = min(indices, key=lambda idx: len(paths[idx]))
            for idx in indices:
                keywords[idx]["isAmbiguousSimpleName"] = True
                keywords[idx]["isShortestForSimpleName"] = idx == shortest_idx
        else:
            keywords[indices[0]]["isAmbiguousSimpleName"] = False
            keywords[indices[0]]["isShortestForSimpleName"] = True

    for i, kw in enumerate(keywords):
        this_path = paths[i]
        for suffix_len in range(1, len(this_path) + 1):
            candidate = ".".join(this_path[-suffix_len:])
            matches = [
                j
                for j, other_path in enumerate(paths)
                if len(other_path) >= suffix_len
                and other_path[-suffix_len:] == this_path[-suffix_len:]
            ]
            if len(matches) == 1 and matches[0] == i:
                kw["minimumPath"] = candidate
                break
        else:
            kw["minimumPath"] = ".".join(this_path)

    return keywords


def argument_names(container: CContainer) -> Dict[str, str]:
    """``{objectPath: name i2run accepts}`` for every parameter in *container*.

    The renderer's whole need, and the reason this module exists: what it emits
    is by construction what i2run's parser was built to accept.
    """
    return {
        kw["path"]: kw["minimumPath"]
        for kw in compute_minimum_paths(leaf_paths(container))
    }
