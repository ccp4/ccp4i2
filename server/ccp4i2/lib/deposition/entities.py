"""Entity bookkeeping for the deposition model, done with gemmi.

The refinement programs write one entity per chain, with no sequence
(``_entity_poly.pdbx_seq_one_letter_code`` is ``?``) and no
``_entity_poly_seq``. The deposited model needs one entity per crystallised
molecule carrying the *full* expressed sequence from the AU contents (tags
and unmodelled loops included -- the PDB records the expressed sequence),
every chain assigned to the sequence it aligns to best, and ``label_seq_id``
numbered on that sequence.

- :func:`rebuild_entities` does that, rewriting only the entity, atom and
  chem_comp categories and keeping ``_refine``, ``_software``,
  ``_struct_ncs_oper`` and the rest exactly as the refinement wrote them.
- :func:`rewrite_entity_poly` replaces ``_entity_poly`` in a file with the
  table gemmi derives from those entities (what the retired PDBe package
  needed after it had appended rows of its own).

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


def residue_kind(polymer_type):
    """(ResidueKind, PolymerType) for an AU-contents polymer type (PROTEIN, RNA, DNA)."""
    return _RESIDUE_KIND.get(str(polymer_type).upper(), _RESIDUE_KIND["PROTEIN"])


def _polymer_entities(sequences):
    entities = []
    for idx, (name, one_letter, polymer_type) in enumerate(sequences, start=1):
        kind, ptype = residue_kind(polymer_type)
        ent = gemmi.Entity(str(idx))
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = ptype
        ent.full_sequence = gemmi.expand_one_letter_sequence(
            "".join(str(one_letter).split()), kind)
        entities.append(ent)
    return entities


def read_model_document(path):
    """A gemmi Document for a model file: the mmCIF as written, or a PDB-format
    file converted (it has no refinement tables to keep)."""
    if str(path).lower().endswith((".cif", ".mmcif", ".cif.gz")):
        return gemmi.cif.read(str(path))
    return gemmi.read_structure(str(path)).make_mmcif_document()


def rebuild_entities(path_in, path_out, sequences):
    """Write ``path_out`` from ``path_in`` with entities derived from ``sequences``.

    ``sequences`` is a list of ``(name, one_letter_sequence, polymer_type)``
    with polymer_type one of PROTEIN, RNA, DNA (the AU-contents vocabulary).
    Returns the gemmi Structure, which :func:`rewrite_entity_poly` needs.
    Raises ValueError if ``sequences`` is empty.
    """
    doc = rebuild_entities_document(path_in, sequences)[0]
    doc.write_file(str(path_out))
    return gemmi.make_structure_from_block(doc[0])


def rebuild_entities_document(path_in, sequences):
    """As :func:`rebuild_entities`, returning ``(document, structure)`` without writing."""
    doc = read_model_document(path_in)
    return doc, rebuild_entities_in(doc, sequences)


MIN_IDENTITY = 0.3


def rebuild_entities_in(doc, sequences, min_identity=MIN_IDENTITY, warnings=None):
    """Rebuild the entities of the first block of ``doc`` in place; returns the Structure.

    A polymer chain is assigned to the AU-contents sequence it aligns to
    best, provided the identity reaches ``min_identity``; a chain below that
    (a peptide ligand, a chain the AU contents forgot) keeps an entity of its
    own with its modelled sequence, and a line is appended to ``warnings``.
    """
    if not sequences:
        raise ValueError("rebuild_entities needs at least one sequence")
    if warnings is None:
        warnings = []
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

    unassigned = []
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        scored = [(identity(ent, poly), ent) for ent in polymers]
        score, best = max(scored, key=lambda t: t[0])
        if score < min_identity:
            # Its own entity, with the sequence as modelled.
            own = gemmi.Entity(str(len(polymers) + len(others) + len(unassigned) + 1))
            own.entity_type = gemmi.EntityType.Polymer
            own.polymer_type = poly.check_polymer_type()
            own.full_sequence = [r.name for r in poly]
            own.subchains = [poly.subchain_id()]
            for res in poly:
                res.entity_id = own.name
            unassigned.append(own)
            warnings.append("chain %s matches no AU-contents sequence (best identity %.0f%%); "
                            "it is deposited as its own entity with its modelled sequence -- "
                            "add it to the AU contents if it is a crystallised molecule"
                            % (chain.name, 100 * score))
            continue
        # Entity.subchains is returned by value from the binding: assign, don't append.
        best.subchains = best.subchains + [poly.subchain_id()]
        for res in poly:
            res.entity_id = best.name
    for ent in others:
        for chain in st[0]:
            for res in chain:
                if res.subchain in ent.subchains:
                    res.entity_id = ent.name

    st.entities = gemmi.EntityList(polymers + unassigned + others)
    st.assign_label_seq_id(force=True)

    groups = gemmi.MmcifOutputGroups(False)
    for name in _ENTITY_GROUPS:
        setattr(groups, name, True)
    st.update_mmcif_block(block, groups)
    status = block.find_mmcif_category("_pdbx_database_status.")
    if status:
        status.erase()
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


def strands_by_entity(block):
    """``{entity_id: [chain ids]}`` from the block's ``_entity_poly`` table."""
    out = {}
    for row in block.find("_entity_poly.", ["entity_id", "pdbx_strand_id"]):
        out[row.str(0)] = [s for s in row.str(1).split(",") if s]
    return out
