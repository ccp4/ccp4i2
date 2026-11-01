"""Sequence database references for the deposited entities.

OneDep's sequence page is populated from ``_struct_ref`` (the reference the
entity's sequence is compared with), ``_struct_ref_seq`` (which stretch of
each chain maps to which stretch of the reference) and
``_struct_ref_seq_dif`` (every residue where construct and reference
disagree, with why). A natural or engineered protein references its UniProt
entry; a designed scaffold with a PDB entry of its own references that; a
de novo design has nothing to reference, and the PDB's convention is then a
``_struct_ref`` to the entry itself (``db_name PDB``) carrying the construct
sequence, so the page is filled rather than empty.

All of it is derived from an alignment of a *segment* of the construct to the
reference sequence, so a chimera (tag + designed domain + linker + natural
domain) gets one reference per segment that has one and the tags fall outside
every alignment.
"""
import gemmi

# Vocabulary of _struct_ref_seq_dif.details used here. The dictionary does not
# enumerate the item, but OneDep's annotators expect these spellings.
DIF_CONFLICT = "conflict"
DIF_ENGINEERED = "engineered mutation"
DIF_INSERTION = "insertion"
DIF_DELETION = "deletion"
DIF_TAG = "expression tag"
DIF_LINKER = "linker"
DIF_CLONING = "cloning artifact"

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def fetch_uniprot(accession, timeout=20):
    """``(db_code, sequence)`` for a UniProt accession, from the REST FASTA.

    Raises ``RuntimeError`` naming the accession when the entry cannot be
    fetched or parsed, so a caller can record the accession and move on.
    """
    import requests  # lazy: execution-time dependency

    acc = str(accession).strip()
    try:
        response = requests.get(UNIPROT_FASTA.format(acc=acc), timeout=timeout)
    except requests.RequestException as err:
        raise RuntimeError("UniProt %s: %s" % (acc, err))
    if response.status_code != 200:
        raise RuntimeError("UniProt %s: HTTP %d" % (acc, response.status_code))
    return parse_fasta(response.text, acc)


RCSB_FASTA = "https://www.rcsb.org/fasta/entry/{pdb_id}/display"


def fetch_pdb_sequence(pdb_id, piece, timeout=20):
    """``(pdb_id, sequence)`` for the polymer entity of a PDB entry that best
    matches ``piece`` (a stretch of the construct), from RCSB's FASTA.

    Raises ``RuntimeError`` when the entry cannot be fetched; the caller then
    records the reference without a sequence.
    """
    import requests  # lazy: execution-time dependency

    code = str(pdb_id).strip().upper()
    try:
        response = requests.get(RCSB_FASTA.format(pdb_id=code), timeout=timeout)
    except requests.RequestException as err:
        raise RuntimeError("PDB %s: %s" % (code, err))
    if response.status_code != 200:
        raise RuntimeError("PDB %s: HTTP %d" % (code, response.status_code))
    best, best_score = None, None
    sequence = ""
    for line in response.text.splitlines() + [">"]:
        if line.startswith(">"):
            if sequence:
                score = sum(1 for p in align_pairs(piece, sequence)
                            if p[0] is not None and p[2] is not None and p[1] == p[3])
                if best_score is None or score > best_score:
                    best, best_score = sequence, score
            sequence = ""
        else:
            sequence += line.strip().upper()
    if not best:
        raise RuntimeError("PDB %s: no polymer sequence in the FASTA" % code)
    return code, best


def parse_fasta(text, accession=""):
    """``(db_code, sequence)`` from a UniProt FASTA (``>sp|ACC|CODE_SPECIES ...``)."""
    lines = [line.strip() for line in str(text).splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError("UniProt %s: not a FASTA response" % accession)
    header = lines[0][1:]
    parts = header.split("|")
    db_code = parts[2].split()[0] if len(parts) >= 3 else header.split()[0]
    sequence = "".join(lines[1:]).upper()
    if not sequence:
        raise RuntimeError("UniProt %s: empty sequence" % accession)
    return db_code, sequence


def align_pairs(construct, reference):
    """Residue-by-residue alignment of two one-letter strings.

    Returns a list of ``(construct_pos, construct_res, ref_pos, ref_res)``
    with 1-based positions; a gap on either side is ``(None, None)`` there.
    End gaps are free, so a reference that covers only the middle of the
    construct does not drag the ends into a bad alignment.
    """
    construct = "".join(str(construct).split()).upper()
    reference = "".join(str(reference).split()).upper()
    scoring = gemmi.AlignmentScoring()
    result = gemmi.align_string_sequences(list(construct), list(reference), [], scoring)
    pairs = []
    ci = ri = 0
    # Expand the CIGAR ("12M3I...") into per-residue operations.
    count = ""
    for ch in result.cigar_str():
        if ch.isdigit():
            count += ch
            continue
        n = int(count or "1")
        count = ""
        for _ in range(n):
            if ch == "M":
                pairs.append((ci + 1, construct[ci], ri + 1, reference[ri]))
                ci += 1
                ri += 1
            elif ch == "I":
                pairs.append((ci + 1, construct[ci], None, None))
                ci += 1
            elif ch == "D":
                pairs.append((None, None, ri + 1, reference[ri]))
                ri += 1
    return pairs


def gapless_segments(pairs):
    """Maximal runs of aligned residues with a constant offset:
    ``[(construct_beg, construct_end, ref_beg, ref_end), ...]``."""
    segments = []
    current = None
    for cpos, cres, rpos, rres in pairs:
        if cpos is None or rpos is None:
            current = None
            continue
        offset = rpos - cpos
        if current and current[3] == offset:
            current[1] = cpos
            continue
        current = [cpos, cpos, rpos, offset]
        segments.append(current)
    return [(c0, c1, r0, r0 + (c1 - c0)) for c0, c1, r0, _ in segments]


def _parse_mutations(text):
    """``{(ref_res, ref_pos, construct_res)}`` from ``"L33H, Q141P"``."""
    out = set()
    for token in str(text or "").replace(";", ",").split(","):
        token = token.strip().upper()
        if len(token) < 3 or not token[1:-1].isdigit():
            continue
        out.add((token[0], int(token[1:-1]), token[-1]))
    return out


def differences(pairs, engineered="", segment_bounds=None):
    """``[(construct_pos, construct_res, ref_pos, ref_res, details), ...]`` for
    every disagreement between construct and reference within the aligned
    region. Construct residues outside the alignment are not reported here:
    they belong to another segment (a tag, a linker, another domain) and are
    described by it.

    ``engineered`` names the intended substitutions (``"L33H"``): a mismatch
    that matches one is an *engineered mutation*, any other is a *conflict*.
    ``segment_bounds`` ``(beg, end)`` limits reporting of construct-only
    residues to insertions inside the segment.
    """
    mutations = _parse_mutations(engineered)
    aligned = [p for p in pairs if p[0] is not None and p[2] is not None]
    if not aligned:
        return []
    first_c, last_c = aligned[0][0], aligned[-1][0]
    out = []
    for cpos, cres, rpos, rres in pairs:
        if cpos is not None and rpos is not None:
            if cres != rres:
                detail = DIF_ENGINEERED if (rres, rpos, cres) in mutations else DIF_CONFLICT
                out.append((cpos, cres, rpos, rres, detail))
        elif cpos is not None:
            if first_c < cpos < last_c:
                out.append((cpos, cres, None, None, DIF_INSERTION))
        else:
            if any(p[2] is not None and p[2] < rpos for p in aligned) and \
                    any(p[2] is not None and p[2] > rpos for p in aligned):
                out.append((None, None, rpos, rres, DIF_DELETION))
    return out


class Reference:
    """One database reference for an entity, over one segment of its construct."""

    def __init__(self, entity_id, db_name, db_code, accession, ref_sequence,
                 construct, beg, end, engineered=""):
        self.entity_id = str(entity_id)
        self.db_name = db_name
        self.db_code = db_code
        self.accession = accession
        self.ref_sequence = "".join(str(ref_sequence).split()).upper()
        self.construct = "".join(str(construct).split()).upper()
        self.beg = int(beg)
        self.end = int(end)
        self.engineered = engineered
        piece = self.construct[self.beg - 1:self.end]
        pairs = align_pairs(piece, self.ref_sequence)
        shift = self.beg - 1
        self.pairs = [(c + shift if c else None, cr, r, rr) for c, cr, r, rr in pairs]
        self.segments = gapless_segments(self.pairs)
        self.diffs = differences(self.pairs, engineered)

    def coverage(self):
        """``(construct_beg, construct_end)`` of the aligned region, or None."""
        if not self.segments:
            return None
        return self.segments[0][0], self.segments[-1][1]


def self_reference(entity_id, construct, entry_id="XXXX"):
    """The PDB convention for a sequence with no database entry (a de novo
    design): a reference to this entry itself, carrying the whole construct."""
    ref = Reference.__new__(Reference)
    ref.entity_id = str(entity_id)
    ref.db_name = "PDB"
    ref.db_code = entry_id
    ref.accession = entry_id
    ref.ref_sequence = "".join(str(construct).split()).upper()
    ref.construct = ref.ref_sequence
    ref.beg, ref.end = 1, len(ref.construct)
    ref.engineered = ""
    ref.pairs = [(i, c, i, c) for i, c in enumerate(ref.construct, start=1)]
    ref.segments = [(1, len(ref.construct), 1, len(ref.construct))] if ref.construct else []
    ref.diffs = []
    return ref


def tag_rows(construct, covered, segments_info):
    """Construct residues no reference covers, as ``_struct_ref_seq_dif``
    details: the segment they fall in says whether they are a tag, a linker or
    a cloning artifact. ``covered`` is a set of construct positions some
    reference aligns; ``segments_info`` is ``[(beg, end, details)]``.

    Returns ``[(construct_pos, construct_res, details)]``.
    """
    out = []
    construct = "".join(str(construct).split()).upper()
    for beg, end, details in segments_info:
        for pos in range(int(beg), int(end) + 1):
            if pos in covered or pos < 1 or pos > len(construct):
                continue
            out.append((pos, construct[pos - 1], details))
    return out


def three_letter(one):
    """Three-letter code for a one-letter residue (``X`` for anything unknown)."""
    try:
        return gemmi.find_tabulated_residue(gemmi.expand_one_letter(one, gemmi.ResidueKind.AA)).name
    except Exception:
        return "UNK"


def reference_tables(references, strands, entry_id="XXXX", extra_difs=()):
    """The three categories as ``{category: (tags, rows)}``.

    ``references`` is a list of :class:`Reference` (any number per entity);
    ``strands`` maps entity id to its chain ids; ``extra_difs`` is
    ``[(entity_id, construct_pos, construct_res, details)]`` for residues no
    reference covers (tags, linkers), reported against the entity's first
    reference so OneDep sees them.
    """
    struct_ref = []
    struct_ref_seq = []
    struct_ref_seq_dif = []
    align_id = 0
    ordinal = 0
    first_align_for = {}
    for ref_id, ref in enumerate(references, start=1):
        align_begin = ref.segments[0][2] if ref.segments else 1
        struct_ref.append((ref_id, ref.db_name, ref.db_code, ref.accession, ref.entity_id,
                           ref.ref_sequence, align_begin))
        for chain in strands.get(ref.entity_id, []):
            for c0, c1, r0, r1 in ref.segments:
                align_id += 1
                first_align_for.setdefault((ref.entity_id, chain), align_id)
                struct_ref_seq.append((align_id, ref_id, entry_id, chain, c0, "?", c1, "?",
                                       ref.accession, r0, "?", r1, "?", c0, c1))
            aid = first_align_for.get((ref.entity_id, chain))
            for cpos, cres, rpos, rres, details in ref.diffs:
                ordinal += 1
                struct_ref_seq_dif.append((
                    aid, entry_id, three_letter(cres) if cres else "?", chain,
                    cpos if cpos else "?", "?", ref.db_name, ref.accession,
                    three_letter(rres) if rres else "?", rpos if rpos else "?",
                    details, cpos if cpos else "?", ordinal))
    for entity_id, cpos, cres, details in extra_difs:
        for chain in strands.get(str(entity_id), []):
            aid = first_align_for.get((str(entity_id), chain))
            if aid is None:
                continue
            ref = next(r for r in references if r.entity_id == str(entity_id))
            ordinal += 1
            struct_ref_seq_dif.append((
                aid, entry_id, three_letter(cres), chain, cpos, "?", ref.db_name,
                ref.accession, "?", "?", details, cpos, ordinal))
    return {
        "_struct_ref": (
            ["id", "db_name", "db_code", "pdbx_db_accession", "entity_id",
             "pdbx_seq_one_letter_code", "pdbx_align_begin"], struct_ref),
        "_struct_ref_seq": (
            ["align_id", "ref_id", "pdbx_PDB_id_code", "pdbx_strand_id", "seq_align_beg",
             "pdbx_seq_align_beg_ins_code", "seq_align_end", "pdbx_seq_align_end_ins_code",
             "pdbx_db_accession", "db_align_beg", "pdbx_db_align_beg_ins_code", "db_align_end",
             "pdbx_db_align_end_ins_code", "pdbx_auth_seq_align_beg", "pdbx_auth_seq_align_end"],
            struct_ref_seq),
        "_struct_ref_seq_dif": (
            ["align_id", "pdbx_pdb_id_code", "mon_id", "pdbx_pdb_strand_id", "seq_num",
             "pdbx_pdb_ins_code", "pdbx_seq_db_name", "pdbx_seq_db_accession_code", "db_mon_id",
             "pdbx_seq_db_seq_num", "details", "pdbx_auth_seq_num", "pdbx_ordinal"],
            struct_ref_seq_dif),
    }
