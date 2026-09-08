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
