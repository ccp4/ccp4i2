"""Server-side harvest helpers shared by the Coot task wrappers.

Unlike api_client (which is loaded into Coot's embedded interpreter and
must stay stdlib-only/py2-py3), this module runs only in the wrapper
under ccp4-python, so it may use gemmi.
"""

from __future__ import absolute_import


def cif_is_restraint_dictionary(path):
    """True if a CIF is a monomer restraint dictionary rather than a
    coordinate model.

    Coot's ligand builder / get-monomer writes a restraint dictionary
    (``_chem_comp_atom`` / ``comp_<code>`` blocks, no ``_atom_site``).
    Harvest code that globs ``*.cif`` would otherwise mis-file it - as a
    coordinate model (coot1), or miss it entirely if it relies on
    filename patterns (coot_rebuild). Told apart by content, so any
    naming the builder uses is handled. Defensive: on any read error,
    returns False (treat as a model / not-a-dict - the prior behaviour).
    """
    try:
        import gemmi

        doc = gemmi.cif.read(str(path))
    except Exception:
        return False

    def has(block, tag):
        try:
            column = block.find_values(tag)
            return column is not None and len(column) > 0
        except Exception:
            return False

    saw_dictionary = False
    for block in doc:
        # A coordinate model is decisive: atom sites mean it is not a dict.
        if has(block, "_atom_site.Cartn_x") or has(block, "_atom_site.id"):
            return False
        if (has(block, "_chem_comp_atom.atom_id")
                or has(block, "_chem_comp_bond.atom_id_1")
                or block.name == "comp_list"
                or block.name.startswith("comp_")):
            saw_dictionary = True
    return saw_dictionary


# ---------------------------------------------------------------------------
# Filing harvested files into output lists (shared by coot1, coot_rebuild
# and the moorhen task)
# ---------------------------------------------------------------------------

import json
import os
from pathlib import Path


def read_drop_metadata(path):
    """Metadata a saver left beside a dropped file, or ``{}``.

    The moorhen task's drop endpoint writes ``output<N>.meta.json`` next to
    ``output<N>.pdb`` carrying the kind and the annotation the user gave.
    Read before filing, because filing may move the file out of the drop
    directory and leave the sidecar behind.
    """
    sidecar = Path(path).with_suffix(".meta.json")
    if not sidecar.is_file():
        return {}
    try:
        with open(sidecar, "r") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def split_models_and_dictionaries(paths):
    """Partition harvested paths into (models, dictionaries) by content.

    A ``.cif`` is a restraint dictionary only if it says so; anything else,
    including a ``.cif`` that is a coordinate model, is a model.
    """
    models = []
    dictionaries = []
    for path in paths:
        path = Path(path)
        if path.suffix == ".cif" and cif_is_restraint_dictionary(path):
            dictionaries.append(path)
        else:
            models.append(path)
    return models, dictionaries


def file_into_list(out_list, paths, work_dir, stem, annotate):
    """File ``paths`` into the ``out_list`` COutputFileList.

    Files outside the work directory are moved to canonical names
    (``<stem>_<i><suffix>``) first. ``annotate(item, path, meta)`` sets the
    metadata; ``meta`` is the sidecar read before the move. Spare slots are
    truncated in place with pop() -- NOT ``out_list.set(slice)``, which
    deep-copies items through CDataFile.get()/set() and drops the
    annotation/subType just set (the gleaner then falls back to the bare
    param name). Returns the number filed.
    """
    work_dir = Path(work_dir)
    index = 0
    for path in paths:
        path = Path(path)
        meta = read_drop_metadata(path)
        if path.parent != work_dir:
            target = work_dir / f"{stem}_{index}{path.suffix}"
            while target.exists():
                target = work_dir / f"{stem}_{index}_{target.stem}{path.suffix}"
            os.replace(path, target)
            path = target
        while index >= len(out_list):
            out_list.append(out_list.makeItem())
        out_list[index].setFullPath(str(path))
        annotate(out_list[index], path, meta)
        index += 1
    while len(out_list) > index:
        out_list.pop()
    return index


def harvest_candidates(work_dir, paths, xyz_list, dict_list,
                       annotate_model, annotate_dict,
                       model_stem="XYZOUT", dict_stem="DICTOUT"):
    """Split ``paths`` into models and dictionaries and file each into its
    output list. Returns ``(n_models, n_dicts)``."""
    models, dictionaries = split_models_and_dictionaries(paths)
    n_models = file_into_list(xyz_list, models, work_dir, model_stem,
                              annotate_model)
    n_dicts = file_into_list(dict_list, dictionaries, work_dir, dict_stem,
                             annotate_dict)
    return n_models, n_dicts


def harvest_drop_directory(work_dir, drop_dir, xyz_list, dict_list,
                           annotate_model, annotate_dict):
    """Harvest the ``output<N>.pdb|cif`` drop contract from ``drop_dir``
    in save order. Returns ``(n_models, n_dicts)``."""
    from . import api_client

    paths = [Path(path) for _number, path
             in api_client.harvestable_outputs(str(drop_dir))]
    return harvest_candidates(work_dir, paths, xyz_list, dict_list,
                              annotate_model, annotate_dict)
