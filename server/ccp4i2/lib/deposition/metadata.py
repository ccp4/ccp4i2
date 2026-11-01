"""The descriptive metadata of a deposition, written into the model mmCIF.

OneDep reads the title, authors, contact, citation, source organisms,
sequence references, data collection, software and assembly from the model
file and pre-fills its form pages from them; whatever is absent is a page the
depositor types by hand. :func:`write_metadata` takes a plain dict (the task
builds it from its parameters; tests build it directly) and writes each
category once, contiguously, through :mod:`mmcif_edit`.

The dict (every key optional unless said):

- ``title`` (required), ``keywords`` (``_struct_keywords.pdbx_keywords``),
  ``keyword_text``;
- ``authors``: ``[{"name": "Noble, M.E.M.", "orcid": "..."}]``;
- ``citation``: ``{"published": bool, "title", "journal", "volume",
  "first_page", "last_page", "year", "doi", "pubmed", "authors": [...]}``
  (unpublished: journal "To Be Published", authors default to the deposition
  authors);
- ``contact``: ``{"salutation", "first_name", "last_name", "role",
  "organization_type", "address_1", "address_2", "city", "state",
  "postal_code", "country", "email", "phone", "orcid"}``;
- ``related``: ``[{"db_name", "db_id", "content_type", "details"}]``;
- ``entities``: one per polymer entity id, ``{"entity_id", "description",
  "common_name", "ec", "mutation", "details", "segments": [...]}``; a segment
  is ``{"beg", "end", "kind", "organism", "taxid", "strain", "gene",
  "host", "host_taxid", "uniprot", "pdb_id", "engineered", "details"}`` with
  ``kind`` one of :data:`SEGMENT_KINDS`;
- ``design_details``: free text on how a designed entity was made;
- ``experiment``: ``{"source", "site", "beamline", "detector", "detector_type",
  "collection_date", "temperature", "wavelength", "monochromator",
  "protocol"}``;
- ``crystal``: ``{"method", "temp", "ph", "details"}``;
- ``software``: ``[{"name", "version", "classification"}]`` added to what the
  refinement program wrote;
- ``phasing``: ``{"method", "starting_model", "details"}``;
- ``assembly``: ``{"mode": "asu"|"author", "oligomeric_count", "operators",
  "details"}``.
"""
import gemmi

from . import struct_ref as sr
from .assembly import assembly_tables, matthews
from .entities import strands_by_entity
from .mmcif_edit import (
    replace_loop,
    replace_pairs,
    set_pairs_in_place,
    split_categories,
    value_of,
)
from .vocab import ASSEMBLY_EVIDENCE, CLASSIFICATIONS, canonical_software_name

# How a segment of a construct came to be. The first-class cases are natural,
# engineered-from-natural, and designed: a world with many de novo proteins
# must not have to invent a UniProt accession or an organism for them.
NATURAL = "natural"              # a natural protein: UniProt reference
ENGINEERED = "engineered"        # variant of a natural protein: UniProt + differences
DESIGNED = "designed"            # de novo design: no sequence database reference
SCAFFOLD = "scaffold"            # a designed scaffold the PDB already holds (a DARPin)
SYNTHETIC = "synthetic_peptide"  # chemically synthesised: _pdbx_entity_src_syn
TAG = "tag"                      # expression / purification tag
LINKER = "linker"
SEGMENT_KINDS = (NATURAL, ENGINEERED, DESIGNED, SCAFFOLD, SYNTHETIC, TAG, LINKER)

SYNTHETIC_CONSTRUCT = ("synthetic construct", 32630)
DEFAULT_HOST = ("Escherichia coli", 562)
DE_NOVO_KEYWORD = "DE NOVO PROTEIN"

GEN_TAGS = ["entity_id", "pdbx_src_id", "pdbx_alt_source_flag", "pdbx_seq_type",
            "pdbx_beg_seq_num", "pdbx_end_seq_num", "pdbx_gene_src_gene",
            "pdbx_gene_src_scientific_name", "pdbx_gene_src_ncbi_taxonomy_id",
            "gene_src_strain", "pdbx_host_org_scientific_name",
            "pdbx_host_org_ncbi_taxonomy_id", "pdbx_description"]
SYN_TAGS = ["entity_id", "pdbx_src_id", "pdbx_alt_source_flag", "pdbx_beg_seq_num",
            "pdbx_end_seq_num", "organism_scientific", "ncbi_taxonomy_id", "details"]

# Every category write_metadata emits; what the dictionary gate checks.
WRITTEN_CATEGORIES = (
    "_struct", "_struct_keywords", "_audit_author", "_citation", "_citation_author",
    "_pdbx_contact_author", "_pdbx_database_related", "_entity", "_entity_name_com",
    "_pdbx_entity_nonpoly", "_entity_src_gen", "_pdbx_entity_src_syn", "_struct_ref",
    "_struct_ref_seq", "_struct_ref_seq_dif", "_diffrn", "_diffrn_source", "_diffrn_detector",
    "_diffrn_radiation", "_diffrn_radiation_wavelength", "_exptl", "_exptl_crystal",
    "_exptl_crystal_grow", "_software", "_refine_hist", "_phasing", "_pdbx_phasing_MR",
    "_pdbx_struct_assembly", "_pdbx_struct_assembly_gen", "_pdbx_struct_oper_list",
    "_pdbx_struct_assembly_auth_evidence", "_pdbx_initial_refinement_model",
    "_pdbx_entity_instance_feature",
)

NONPOLYMER_NAMES = {
    "HOH": "water", "CL": "CHLORIDE ION", "NA": "SODIUM ION", "K": "POTASSIUM ION",
    "MG": "MAGNESIUM ION", "CA": "CALCIUM ION", "ZN": "ZINC ION", "MN": "MANGANESE (II) ION",
    "FE": "FE (III) ION", "NI": "NICKEL (II) ION", "CU": "COPPER (II) ION", "SO4": "SULFATE ION",
    "PO4": "PHOSPHATE ION", "GOL": "GLYCEROL", "EDO": "1,2-ETHANEDIOL", "PEG": "DI(HYDROXYETHYL)ETHER",
    "ACT": "ACETATE ION", "DMS": "DIMETHYL SULFOXIDE", "NAG": "2-acetamido-2-deoxy-beta-D-glucopyranose",
    "ATP": "ADENOSINE-5'-TRIPHOSPHATE", "ADP": "ADENOSINE-5'-DIPHOSPHATE", "GDP": "GUANOSINE-5'-DIPHOSPHATE",
    "GTP": "GUANOSINE-5'-TRIPHOSPHATE", "NAD": "NICOTINAMIDE-ADENINE-DINUCLEOTIDE",
    "NAP": "NADP NICOTINAMIDE-ADENINE-DINUCLEOTIDE PHOSPHATE", "FAD": "FLAVIN-ADENINE DINUCLEOTIDE",
    "HEM": "PROTOPORPHYRIN IX CONTAINING FE", "TRS": "2-AMINO-2-HYDROXYMETHYL-PROPANE-1,3-DIOL",
    "MPD": "(4S)-2-METHYL-2,4-PENTANEDIOL", "IMD": "IMIDAZOLE", "BME": "BETA-MERCAPTOETHANOL",
}


def _get(d, key, default=None):
    value = (d or {}).get(key, default)
    if value is None:
        return default
    if isinstance(value, str) and not value.strip():
        return default
    return value


def _one_letter(entity):
    return gemmi.one_letter_code([gemmi.Entity.first_mon(m) for m in entity.full_sequence]).upper()


def entity_sequences(structure):
    """``{entity_id: one-letter sequence}`` for the polymer entities."""
    return {ent.name: _one_letter(ent) for ent in structure.entities
            if ent.entity_type == gemmi.EntityType.Polymer}


def _segments_for(entity_spec, length):
    """The entity's segments, defaulting to one segment over the whole
    sequence; each with integer bounds inside the sequence."""
    segments = []
    for seg in (entity_spec or {}).get("segments") or []:
        beg = int(_get(seg, "beg", 1))
        end = int(_get(seg, "end", length))
        beg, end = max(1, beg), min(length, end)
        if end < beg:
            continue
        kind = str(_get(seg, "kind", NATURAL)).strip().lower()
        if kind not in SEGMENT_KINDS:
            kind = NATURAL
        copy = dict(seg)
        copy.update({"beg": beg, "end": end, "kind": kind})
        segments.append(copy)
    if not segments:
        segments = [{"beg": 1, "end": length, "kind": NATURAL}]
    return segments


def _seq_type(seg, length):
    if seg["kind"] == TAG:
        if seg["beg"] == 1:
            return "N-terminal tag"
        if seg["end"] == length:
            return "C-terminal tag"
        return "Linker"
    if seg["kind"] == LINKER:
        return "Linker"
    return "Biological sequence"


def _src_method(segments):
    kinds = {s["kind"] for s in segments}
    if SYNTHETIC in kinds:
        return "syn"
    if kinds <= {NATURAL}:
        return "man"
    return "man"


def _source_rows(entity_id, segments, length, src_method):
    gen, syn = [], []
    for n, seg in enumerate(segments, start=1):
        kind = seg["kind"]
        if src_method == "syn":
            syn.append((entity_id, n, "sample", seg["beg"], seg["end"],
                        _get(seg, "organism", SYNTHETIC_CONSTRUCT[0]),
                        _get(seg, "taxid", SYNTHETIC_CONSTRUCT[1]), _get(seg, "details")))
            continue
        if kind in (NATURAL, ENGINEERED):
            organism, taxid = _get(seg, "organism"), _get(seg, "taxid")
        else:
            organism, taxid = SYNTHETIC_CONSTRUCT
        gen.append((entity_id, n, "sample", _seq_type(seg, length), seg["beg"], seg["end"],
                    _get(seg, "gene") if kind in (NATURAL, ENGINEERED) else None,
                    organism, taxid,
                    _get(seg, "strain") if kind in (NATURAL, ENGINEERED) else None,
                    _get(seg, "host", DEFAULT_HOST[0]), _get(seg, "host_taxid", DEFAULT_HOST[1]),
                    _get(seg, "details")))
    return gen, syn


def _references_for(entity_id, sequence, segments, fetch, warnings, entry_id,
                    fetch_pdb=sr.fetch_pdb_sequence):
    """``(references, extra_difs, covered)`` for one entity."""
    refs, extra = [], []
    covered = set()
    length = len(sequence)
    for seg in segments:
        kind = seg["kind"]
        piece_detail = None
        if kind in (NATURAL, ENGINEERED) and _get(seg, "uniprot"):
            acc = str(seg["uniprot"]).strip()
            try:
                db_code, ref_seq = fetch(acc)
            except Exception as err:  # network, unknown accession
                warnings.append("%s; the reference is recorded without its sequence" % err)
                db_code, ref_seq = acc, None
            if ref_seq:
                ref = sr.Reference(entity_id, "UNP", db_code, acc, ref_seq, sequence,
                                   seg["beg"], seg["end"], _get(seg, "engineered", ""))
                refs.append(ref)
                for c0, c1, _, _ in ref.segments:
                    covered.update(range(c0, c1 + 1))
            else:
                refs.append(_bare_reference(entity_id, "UNP", db_code, acc, sequence, seg))
                covered.update(range(seg["beg"], seg["end"] + 1))
        elif kind == SCAFFOLD and _get(seg, "pdb_id"):
            pdb_id = str(seg["pdb_id"]).strip().upper()
            piece = sequence[seg["beg"] - 1:seg["end"]]
            try:
                _, ref_seq = fetch_pdb(pdb_id, piece)
            except Exception as err:
                warnings.append("%s; the scaffold reference is recorded without its sequence" % err)
                ref_seq = None
            if ref_seq:
                ref = sr.Reference(entity_id, "PDB", pdb_id, pdb_id, ref_seq, sequence,
                                   seg["beg"], seg["end"], _get(seg, "engineered", ""))
                refs.append(ref)
                for c0, c1, _, _ in ref.segments:
                    covered.update(range(c0, c1 + 1))
            else:
                refs.append(_bare_reference(entity_id, "PDB", pdb_id, pdb_id, sequence, seg))
                covered.update(range(seg["beg"], seg["end"] + 1))
        elif kind == TAG:
            piece_detail = sr.DIF_TAG
        elif kind == LINKER:
            piece_detail = sr.DIF_LINKER
        if piece_detail:
            extra.append((seg["beg"], seg["end"], piece_detail))
    if not refs:
        # A de novo design (or an entity nobody gave a reference for): the PDB
        # convention is a reference to the entry itself with the whole sequence.
        refs.append(sr.self_reference(entity_id, sequence, entry_id))
        return refs, [], set(range(1, length + 1))
    difs = [(entity_id, pos, res, det) for pos, res, det in sr.tag_rows(sequence, covered, extra)]
    return refs, difs, covered


def _bare_reference(entity_id, db_name, db_code, accession, sequence, seg):
    """A reference whose sequence could not be fetched: the segment maps to
    the database entry one-to-one from its start, with no differences."""
    ref = sr.Reference.__new__(sr.Reference)
    ref.entity_id = str(entity_id)
    ref.db_name, ref.db_code, ref.accession = db_name, db_code, accession
    piece = sequence[seg["beg"] - 1:seg["end"]]
    ref.ref_sequence = piece
    ref.construct = sequence
    ref.beg, ref.end = seg["beg"], seg["end"]
    ref.engineered = ""
    db_beg = int(_get(seg, "db_beg", 1))
    ref.pairs = [(seg["beg"] + i, c, db_beg + i, c) for i, c in enumerate(piece)]
    ref.segments = [(seg["beg"], seg["end"], db_beg, db_beg + len(piece) - 1)] if piece else []
    ref.diffs = []
    return ref


def default_keywords(entity_specs):
    """``DE NOVO PROTEIN`` when every molecule is designed (tags and linkers
    aside); nothing otherwise, because the classification is one term from
    OneDep's list and a chimera's function class is the depositor's call."""
    kinds = set()
    for ent in entity_specs or []:
        for seg in ent.get("segments") or []:
            kind = str(seg.get("kind", NATURAL)).lower()
            if kind not in (TAG, LINKER):
                kinds.add(kind)
    if kinds and kinds <= {DESIGNED, SCAFFOLD}:
        return DE_NOVO_KEYWORD
    return None


def write_metadata(block, structure, spec, entry_id="XXXX", fetch=sr.fetch_uniprot,
                   fetch_pdb=sr.fetch_pdb_sequence):
    """Write the descriptive categories into ``block`` for ``structure``.

    Returns a list of warnings (things recorded incompletely). Raises
    ``ValueError`` when the result would have a split category, which OneDep
    rejects outright.
    """
    warnings = []
    spec = spec or {}
    entry_id = str(_get(spec, "entry_id", entry_id))
    sequences = entity_sequences(structure)
    strands = strands_by_entity(block)
    entity_specs = {str(_get(e, "entity_id", i + 1)): e
                    for i, e in enumerate(spec.get("entities") or [])}

    # --- title, keywords ---------------------------------------------------
    title = _get(spec, "title")
    if not title:
        raise ValueError("a deposition needs a title")
    replace_pairs(block, "_struct", [("entry_id", entry_id), ("title", title),
                                     ("pdbx_CASP_flag", "N")])
    keywords = _get(spec, "keywords") or default_keywords(spec.get("entities"))
    if not keywords:
        raise ValueError("a classification (one term from OneDep's list, such as "
                         "OXIDOREDUCTASE or DE NOVO PROTEIN) is needed")
    keywords = str(keywords).strip().upper()
    if keywords not in CLASSIFICATIONS:
        raise ValueError("classification %r is not one term of OneDep's list; put the rest "
                         "in the keywords text" % keywords)
    replace_pairs(block, "_struct_keywords", [("entry_id", entry_id), ("pdbx_keywords", keywords),
                                              ("text", _get(spec, "keyword_text"))])

    # --- people -----------------------------------------------------------
    authors = [a for a in (spec.get("authors") or []) if _get(a, "name")]
    if authors:
        replace_loop(block, "_audit_author", ["name", "pdbx_ordinal", "identifier_ORCID"],
                     [(a["name"], i, _get(a, "orcid")) for i, a in enumerate(authors, start=1)])
    citation = spec.get("citation") or {}
    published = bool(_get(citation, "published", False))
    replace_pairs(block, "_citation", [
        ("id", "primary"),
        ("title", _get(citation, "title", title)),
        ("journal_abbrev", _get(citation, "journal", "To Be Published") if published else "To Be Published"),
        ("journal_volume", _get(citation, "volume") if published else None),
        ("page_first", _get(citation, "first_page") if published else None),
        ("page_last", _get(citation, "last_page") if published else None),
        ("year", _get(citation, "year") if published else None),
        ("pdbx_database_id_DOI", _get(citation, "doi") if published else None),
        ("pdbx_database_id_PubMed", _get(citation, "pubmed") if published else None),
        ("unpublished_flag", "N" if published else "Y"),
    ])
    cit_authors = [a for a in (citation.get("authors") or []) if _get(a, "name")] or authors
    if cit_authors:
        replace_loop(block, "_citation_author", ["citation_id", "name", "ordinal", "identifier_ORCID"],
                     [("primary", a["name"], i, _get(a, "orcid")) for i, a in enumerate(cit_authors, 1)])
    contact = spec.get("contact") or {}
    if _get(contact, "last_name"):
        replace_loop(block, "_pdbx_contact_author", [
            "id", "name_salutation", "name_first", "name_last", "role", "organization_type",
            "address_1", "address_2", "city", "state_province", "postal_code", "country",
            "email", "phone", "identifier_ORCID"], [(
                1, _get(contact, "salutation"), _get(contact, "first_name"), contact["last_name"],
                _get(contact, "role", "principal investigator/group leader"),
                _get(contact, "organization_type", "academic"),
                _get(contact, "address_1"), _get(contact, "address_2"), _get(contact, "city"),
                _get(contact, "state"), _get(contact, "postal_code"), _get(contact, "country"),
                _get(contact, "email"), _get(contact, "phone"), _get(contact, "orcid"))])
    related = [r for r in (spec.get("related") or []) if _get(r, "db_id")]
    if related:
        replace_loop(block, "_pdbx_database_related", ["db_name", "db_id", "content_type", "details"],
                     [(_get(r, "db_name", "PDB"), str(r["db_id"]).upper(),
                       _get(r, "content_type", "unspecified"), _get(r, "details")) for r in related])

    # --- entities and their sources -------------------------------------
    entity_rows, nonpoly_rows, names_rows = [], [], []
    gen_rows, syn_rows = [], []
    references, extra_difs = [], []
    for ent in structure.entities:
        eid = ent.name
        copies = max(1, len(ent.subchains))
        if ent.entity_type == gemmi.EntityType.Polymer:
            espec = entity_specs.get(eid, {})
            sequence = sequences[eid]
            segments = _segments_for(espec, len(sequence))
            method = _src_method(segments)
            description = _get(espec, "description", "polymer %s" % eid)
            entity_rows.append((eid, "polymer", method, description, copies,
                                _get(espec, "ec"), _get(espec, "mutation"), None,
                                _get(espec, "details", _get(spec, "design_details") if
                                     any(s["kind"] in (DESIGNED, SCAFFOLD) for s in segments) else None)))
            if _get(espec, "common_name"):
                names_rows.append((eid, espec["common_name"]))
            gen, syn = _source_rows(eid, segments, len(sequence), method)
            gen_rows.extend(gen)
            syn_rows.extend(syn)
            refs, difs, _ = _references_for(eid, sequence, segments, fetch, warnings, entry_id,
                                            fetch_pdb)
            references.extend(refs)
            extra_difs.extend(difs)
        else:
            comp = ent.subchains and _comp_of(structure, ent) or "?"
            kind = "water" if ent.entity_type == gemmi.EntityType.Water else "non-polymer"
            name = NONPOLYMER_NAMES.get(comp, comp)
            entity_rows.append((eid, kind, "syn" if kind != "water" else "nat", name, copies,
                                None, None, None, None))
            if kind == "non-polymer":
                nonpoly_rows.append((eid, name, comp))
    replace_loop(block, "_entity", ["id", "type", "src_method", "pdbx_description",
                                    "pdbx_number_of_molecules", "pdbx_ec", "pdbx_mutation",
                                    "pdbx_fragment", "details"], entity_rows)
    replace_loop(block, "_entity_name_com", ["entity_id", "name"], names_rows)
    replace_loop(block, "_pdbx_entity_nonpoly", ["entity_id", "name", "comp_id"], nonpoly_rows)
    replace_loop(block, "_entity_src_gen", GEN_TAGS, gen_rows)
    replace_loop(block, "_pdbx_entity_src_syn", SYN_TAGS, syn_rows)
    tables = sr.reference_tables(references, strands, entry_id, extra_difs,
                                 auth_map=sr.auth_numbering(structure))
    for category, (tags, rows) in tables.items():
        replace_loop(block, category, tags, rows)

    # --- experiment -----------------------------------------------------
    exp = spec.get("experiment") or {}
    wavelength_id = value_of(block, "_diffrn_radiation_wavelength.id")
    wavelength = _get(exp, "wavelength") or value_of(block, "_diffrn_radiation_wavelength.wavelength")
    n_wavelengths = len(block.find_mmcif_category("_diffrn_radiation_wavelength."))
    if n_wavelengths > 1:
        # MAD: every wavelength stays as the refinement program wrote it.
        if _get(exp, "wavelength"):
            warnings.append("the model carries %d wavelengths; the form's %s was not written over them"
                            % (n_wavelengths, _get(exp, "wavelength")))
    elif wavelength_id is None or _get(exp, "wavelength"):
        replace_pairs(block, "_diffrn_radiation_wavelength", [("id", wavelength_id or 1),
                                                              ("wavelength", wavelength)])
        wavelength_id = wavelength_id or "1"
    if not wavelength:
        warnings.append("no wavelength: the model has none and the form gave none")
    replace_pairs(block, "_diffrn", [("id", 1), ("crystal_id", 1),
                                     ("ambient_temp", _get(exp, "temperature")),
                                     ("ambient_temp_details", None)])
    source = _get(exp, "source", "SYNCHROTRON")
    site, beamline = _get(exp, "site"), _get(exp, "beamline")
    source_type = _get(exp, "source_type") or (
        "%s BEAMLINE %s" % (site.upper(), beamline) if site and beamline and source == "SYNCHROTRON"
        else None)
    replace_pairs(block, "_diffrn_source", [
        ("diffrn_id", 1), ("source", source), ("type", source_type),
        ("pdbx_synchrotron_site", site if source == "SYNCHROTRON" else None),
        ("pdbx_synchrotron_beamline", beamline if source == "SYNCHROTRON" else None),
        ("pdbx_wavelength_list", wavelength)])
    replace_pairs(block, "_diffrn_detector", [
        ("diffrn_id", 1), ("detector", _get(exp, "detector")), ("type", _get(exp, "detector_type")),
        ("pdbx_collection_date", _get(exp, "collection_date")), ("details", None)])
    replace_pairs(block, "_diffrn_radiation", [
        ("diffrn_id", 1), ("wavelength_id", wavelength_id),
        ("pdbx_monochromatic_or_laue_m_l", "M"), ("monochromator", _get(exp, "monochromator")),
        ("pdbx_diffrn_protocol", _get(exp, "protocol", "SINGLE WAVELENGTH")),
        ("pdbx_scattering_type", "x-ray")])
    vm, solvent = matthews(structure)
    replace_pairs(block, "_exptl_crystal", [("id", 1), ("density_Matthews", vm),
                                            ("density_percent_sol", solvent), ("description", None)])
    crystal = spec.get("crystal") or {}
    replace_pairs(block, "_exptl_crystal_grow", [
        ("crystal_id", 1), ("method", _get(crystal, "method")), ("temp", _get(crystal, "temp")),
        ("pH", _get(crystal, "ph")), ("pdbx_details", _get(crystal, "details"))])
    set_pairs_in_place(block, "_exptl", [("entry_id", entry_id), ("method", "X-RAY DIFFRACTION"),
                                         ("crystals_number", 1)])

    # --- software -------------------------------------------------------
    existing = []
    for row in block.find("_software.", ["name", "version", "classification"]):
        existing.append({"name": row.str(0), "version": row.str(1), "classification": row.str(2)})
    rows, seen = [], set()
    for sw in list(spec.get("software") or []) + existing:
        name = _get(sw, "name")
        if not name:
            continue
        canonical = canonical_software_name(name)
        if canonical is None:
            # OneDep drops a name outside its software list silently; say so instead.
            warnings.append("software %r is not in OneDep's software list and was not "
                            "written" % name)
            continue
        key = canonical.lower()   # the form's row wins over the file's (it comes first)
        if key in seen:
            continue
        seen.add(key)
        rows.append((len(rows) + 1, canonical, _get(sw, "version"), _get(sw, "classification", "other")))
    if rows:
        replace_loop(block, "_software", ["pdbx_ordinal", "name", "version", "classification"], rows)

    # --- refinement, phasing --------------------------------------------
    n_atoms = structure[0].count_atom_sites()
    n_water = sum(1 for ch in structure[0] for r in ch if r.is_water() for _ in r)
    n_polymer = sum(1 for ch in structure[0] for r in ch.get_polymer() for _ in r)
    n_ligand = n_atoms - n_water - n_polymer
    replace_pairs(block, "_refine_hist", [
        ("pdbx_refine_id", "X-RAY DIFFRACTION"), ("cycle_id", "LAST"),
        ("d_res_high", gemmi.cif.as_string(block.find_value("_refine.ls_d_res_high") or "?")),
        ("d_res_low", gemmi.cif.as_string(block.find_value("_refine.ls_d_res_low") or "?")),
        ("pdbx_number_atoms_protein", n_polymer),
        ("pdbx_number_atoms_ligand", n_ligand), ("pdbx_number_atoms_solvent", n_water),
        ("number_atoms_total", n_atoms)])
    phasing = spec.get("phasing") or {}
    method = _get(phasing, "method")
    models = [m for m in (phasing.get("starting_models") or []) if _get(m, "accession") or _get(m, "details")]
    polymer_ids = [e.name for e in structure.entities if e.entity_type == gemmi.EntityType.Polymer]
    model_rows = []
    for i, m in enumerate(models, start=1):
        model_rows.append((i, _get(m, "entity_ids") or ",".join(polymer_ids),
                           _get(m, "type", "experimental model"), _get(m, "source", "PDB"),
                           _get(m, "accession"), _get(m, "details")))
    replace_loop(block, "_pdbx_initial_refinement_model",
                 ["id", "entity_id_list", "type", "source_name", "accession_code", "details"], model_rows)
    summary = _get(phasing, "starting_model") or "; ".join(
        _get(m, "accession") or _get(m, "details") for m in models) or None
    refine_items = [("pdbx_starting_model", summary)]
    if method:
        replace_pairs(block, "_phasing", [("method", method)])
        refine_items.append(("pdbx_method_to_determine_struct",
                             "MOLECULAR REPLACEMENT" if method == "MR" else method))
        if method == "MR":
            replace_pairs(block, "_pdbx_phasing_MR", [
                ("entry_id", entry_id), ("method_rotation", None), ("method_translation", None),
                ("model_details", _get(phasing, "details", summary))])
    if block.find_mmcif_category("_refine."):
        set_pairs_in_place(block, "_refine", [(k, v) for k, v in refine_items if v])

    # --- assembly -------------------------------------------------------
    asym_ids = [row.str(0) for row in block.find("_struct_asym.", ["id"])]
    asm = spec.get("assembly") or {}
    # Triplets contain commas, so several operators are separated by ";"
    # (or given as a list).
    raw = _get(asm, "operators", "")
    operators = [t.strip() for t in (raw if isinstance(raw, (list, tuple)) else str(raw).split(";"))
                 if str(t).strip()]
    tables = assembly_tables(structure, asym_ids, _get(asm, "mode", "asu"),
                             _get(asm, "oligomeric_count"), operators)
    for category, (tags, rows) in tables.items():
        replace_loop(block, category, tags, rows)
    evidence = str(_get(asm, "evidence", "none")).strip()
    if evidence not in ASSEMBLY_EVIDENCE:
        warnings.append("assembly evidence %r is not in the dictionary's list; wrote 'none'"
                        % evidence)
        evidence = "none"
    replace_loop(block, "_pdbx_struct_assembly_auth_evidence",
                 ["id", "assembly_id", "experimental_support", "details"],
                 [(1, 1, evidence, _get(asm, "details"))])

    # --- ligand of interest -----------------------------------------------
    ligand = str(_get(spec, "ligand_of_interest", "") or "").strip().upper()
    feature_rows = []
    if ligand:
        for chain in structure[0]:
            for res in chain:
                if res.name.upper() == ligand:
                    feature_rows.append((len(feature_rows) + 1, "SUBJECT OF INVESTIGATION",
                                         res.subchain, chain.name, res.seqid.num,
                                         res.label_seq if res.label_seq is not None else None,
                                         res.name, res.name, None))
        if not feature_rows:
            warnings.append("ligand of interest %s is not in the model; nothing written" % ligand)
    replace_loop(block, "_pdbx_entity_instance_feature",
                 ["ordinal", "feature_type", "asym_id", "auth_asym_id", "auth_seq_num", "seq_num",
                  "comp_id", "auth_comp_id", "details"], feature_rows)

    split = split_categories(block)
    if split:
        raise ValueError("categories written in two places (OneDep rejects this): %s"
                         % ", ".join(split))
    return warnings


def _comp_of(structure, entity):
    for chain in structure[0]:
        for res in chain:
            if res.subchain in entity.subchains:
                return res.name
    return "?"
