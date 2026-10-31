"""Entity bookkeeping for the deposition model, done with gemmi.

The refinement programs write one entity per chain, with no sequence
(``_entity_poly.pdbx_seq_one_letter_code`` is ``?``) and no
``_entity_poly_seq``. The PDBe ``adding_stats_to_mmcif`` package then adds
``_entity_poly`` rows *by appending* to whatever loop is there, one entity at a
time, re-reading the loop between appends. With a single entity that gives a
stray duplicate row; with six (a hexamer refined by servalcat) it gave a
23,000-line ``_entity_poly`` loop with the sequence in the wrong column, which
OneDep cannot use.

So the task does the entity work itself, before and after the package runs:

- :func:`rebuild_entities` reads the refinement's mmCIF, makes one polymer
  entity per sequence in the AU contents (the *full* construct, tags and
  unmodelled loops included -- the PDB records the expressed sequence),
  assigns every polymer chain to the sequence it aligns to best, re-derives
  ``label_seq_id`` from that sequence, and rewrites only the entity, atom and
  chem_comp categories, keeping ``_refine``, ``_software``,
  ``_struct_ncs_oper`` and the rest exactly as the refinement wrote them.
- :func:`rewrite_entity_poly` replaces ``_entity_poly`` in the package's
  output with the table gemmi derives from those entities.

The ``_pdbx_database_status`` category is dropped: refinement programs copy it
from whatever PDB file the model started life as, so it carries a deposition
date that is not this entry's.
"""
import gemmi

_RESIDUE_KIND = {
    "PROTEIN": (gemmi.ResidueKind.AA, gemmi.PolymerType.PeptideL),
    "RNA": (gemmi.ResidueKind.RNA, gemmi.PolymerType.Rna),
    "DNA": (gemmi.ResidueKind.DNA, gemmi.PolymerType.Dna),
}

# Output groups that rebuild_entities rewrites; everything else in the block is left alone.
_ENTITY_GROUPS = ("atoms", "entity", "entity_poly", "entity_poly_seq", "struct_asym", "chem_comp")


def _polymer_entities(sequences):
    entities = []
    for idx, (name, one_letter, polymer_type) in enumerate(sequences, start=1):
        kind, ptype = _RESIDUE_KIND.get(str(polymer_type).upper(), _RESIDUE_KIND["PROTEIN"])
        ent = gemmi.Entity(str(idx))
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = ptype
        ent.full_sequence = gemmi.expand_one_letter_sequence(
            "".join(str(one_letter).split()), kind)
        entities.append(ent)
    return entities


def rebuild_entities(path_in, path_out, sequences):
    """Write ``path_out`` from ``path_in`` with entities derived from ``sequences``.

    ``sequences`` is a list of ``(name, one_letter_sequence, polymer_type)``
    with polymer_type one of PROTEIN, RNA, DNA (the AU-contents vocabulary).
    Returns the gemmi Structure, which :func:`rewrite_entity_poly` needs.
    Raises ValueError if ``sequences`` is empty.
    """
    if not sequences:
        raise ValueError("rebuild_entities needs at least one sequence")
    if str(path_in).lower().endswith((".cif", ".mmcif", ".cif.gz")):
        doc = gemmi.cif.read(str(path_in))
    else:
        # A PDB-format model has no refinement tables to keep; just convert it.
        doc = gemmi.read_structure(str(path_in)).make_mmcif_document()
    block = doc[0]
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    st.assign_subchains(force=True)

    polymers = _polymer_entities(sequences)
    others = [e for e in st.entities if e.entity_type != gemmi.EntityType.Polymer]
    for i, ent in enumerate(others, start=len(polymers) + 1):
        ent.name = str(i)

    def identity(ent, poly):
        alignment = gemmi.align_sequence_to_polymer(ent.full_sequence, poly, ent.polymer_type)
        return alignment.calculate_identity()

    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        best = max(polymers, key=lambda ent: identity(ent, poly))
        # Entity.subchains is returned by value from the binding: assign, don't append.
        best.subchains = best.subchains + [poly.subchain_id()]
        for res in poly:
            res.entity_id = best.name
    for ent in others:
        for chain in st[0]:
            for res in chain:
                if res.subchain in ent.subchains:
                    res.entity_id = ent.name

    st.entities = gemmi.EntityList(polymers + others)
    st.assign_label_seq_id(force=True)

    groups = gemmi.MmcifOutputGroups(False)
    for name in _ENTITY_GROUPS:
        setattr(groups, name, True)
    st.update_mmcif_block(block, groups)
    status = block.find_mmcif_category("_pdbx_database_status.")
    if status:
        status.erase()
    doc.write_file(str(path_out))
    return st


def rewrite_entity_poly(path, structure):
    """Replace ``_entity_poly`` in the mmCIF at ``path`` with gemmi's table for ``structure``."""
    doc = gemmi.cif.read(str(path))
    block = doc[0]
    groups = gemmi.MmcifOutputGroups(False)
    groups.entity_poly = True
    structure.update_mmcif_block(block, groups)
    doc.write_file(str(path))


def entity_summary(path):
    """From an mmCIF: (number of polymer entities, {entity_id: [strand ids]},
    {entity_id: sequence length})."""
    block = gemmi.cif.read(str(path))[0]
    strands, lengths = {}, {}
    tags = ["entity_id", "pdbx_strand_id", "pdbx_seq_one_letter_code"]
    for row in block.find("_entity_poly.", tags):
        strands[row.str(0)] = [s for s in row.str(1).split(",") if s]
        lengths[row.str(0)] = len("".join(row.str(2).split()))
    return len(strands), strands, lengths
