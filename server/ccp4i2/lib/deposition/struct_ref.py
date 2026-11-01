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


MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -1, -4, -1


def align_pairs(construct, reference):
    """Residue-by-residue alignment of two one-letter strings.

    Returns a list of ``(construct_pos, construct_res, ref_pos, ref_res)``
    with 1-based positions; a gap on either side is ``(None, None)`` there.

    Gotoh (affine gaps) with free end gaps on both sequences: a reference
    that covers only the middle of the construct is not dragged over the
    flanks, a construct that stops before the reference does is not charged
    for what it lacks, and one long insertion beats a scatter of short ones
    (so a StrepTag inserted in place of two residues reads as fifteen
    inserted residues plus two substitutions, as the archive records it,
    rather than as single matching residues fished out of the tag).
    """
    a = "".join(str(construct).split()).upper()
    b = "".join(str(reference).split()).upper()
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return [(i + 1, c, None, None) for i, c in enumerate(a)] + \
               [(None, None, j + 1, r) for j, r in enumerate(b)]
    neg = float("-inf")
    # M: a[i] aligned to b[j]; X: gap in b (construct-only); Y: gap in a (reference-only)
    M = [[neg] * (m + 1) for _ in range(n + 1)]
    X = [[neg] * (m + 1) for _ in range(n + 1)]
    Y = [[neg] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0.0
    for i in range(1, n + 1):
        X[i][0] = 0.0   # free leading gap in the reference
    for j in range(1, m + 1):
        Y[0][j] = 0.0   # free leading gap in the construct
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sub = MATCH if a[i - 1] == b[j - 1] else MISMATCH
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + sub
            # trailing gaps are free: at the last column/row an extension costs nothing
            xo, xe = (0.0, 0.0) if j == m else (GAP_OPEN, GAP_EXTEND)
            yo, ye = (0.0, 0.0) if i == n else (GAP_OPEN, GAP_EXTEND)
            X[i][j] = max(M[i - 1][j] + xo, X[i - 1][j] + xe, Y[i - 1][j] + xo)
            Y[i][j] = max(M[i][j - 1] + yo, Y[i][j - 1] + ye, X[i][j - 1] + yo)
    i, j = n, m
    state = max(("M", M[n][m]), ("X", X[n][m]), ("Y", Y[n][m]), key=lambda t: t[1])[0]
    pairs = []
    while i > 0 or j > 0:
        if i == 0:
            pairs.append((None, None, j, b[j - 1]))
            j -= 1
            continue
        if j == 0:
            pairs.append((i, a[i - 1], None, None))
            i -= 1
            continue
        if state == "M":
            pairs.append((i, a[i - 1], j, b[j - 1]))
            sub = MATCH if a[i - 1] == b[j - 1] else MISMATCH
            target = M[i][j] - sub
            i, j = i - 1, j - 1
            state = "M" if M[i][j] == target else ("X" if X[i][j] == target else "Y")
        elif state == "X":
            pairs.append((i, a[i - 1], None, None))
            xo, xe = (0.0, 0.0) if j == m else (GAP_OPEN, GAP_EXTEND)
            target = X[i][j]
            i -= 1
            if M[i][j] + xo == target:
                state = "M"
            elif X[i][j] + xe == target:
                state = "X"
            else:
                state = "Y"
        else:
            pairs.append((None, None, j, b[j - 1]))
            yo, ye = (0.0, 0.0) if i == n else (GAP_OPEN, GAP_EXTEND)
            target = Y[i][j]
            j -= 1
            if M[i][j] + yo == target:
                state = "M"
            elif Y[i][j] + ye == target:
                state = "Y"
            else:
                state = "X"
    pairs.reverse()
    return pairs


MIN_FLANK = 3


def trim_flanks(pairs, min_flank=MIN_FLANK):
    """The same alignment with a first or last aligned run shorter than
    ``min_flank`` residues unaligned: a lone matching residue at the edge of a
    gap is a coincidence of the alignment, not a mapping, and its residues
    are better reported as a cloning artifact than as a one-residue segment."""
    pairs = list(pairs)
    for side in (1, -1):
        seq = pairs if side == 1 else pairs[::-1]
        run, i = 0, 0
        while i < len(seq) and (seq[i][0] is None or seq[i][2] is None):
            i += 1
        start = i
        while i < len(seq) and seq[i][0] is not None and seq[i][2] is not None:
            run += 1
            i += 1
        # the run must be followed by a gap for it to count as a flank piece
        if 0 < run < min_flank and i < len(seq):
            fixed = []
            for k, p in enumerate(seq):
                if start <= k < start + run:
                    fixed.append((p[0], p[1], None, None))
                    fixed.append((None, None, p[2], p[3]))
                else:
                    fixed.append(p)
            seq = fixed
        pairs = seq if side == 1 else seq[::-1]
    return unalign_terminal_mismatches(pairs)


def unalign_terminal_mismatches(pairs):
    """Mismatched residues before the first identical pair and after the
    last are not a mapping to the reference: the construct's own junction
    (LQ where the gene reads MK; the RS a vector leaves behind). They come
    back as construct-only and reference-only positions."""
    aligned = [k for k, p in enumerate(pairs) if p[0] is not None and p[2] is not None]
    matches = [k for k in aligned if pairs[k][1] == pairs[k][3]]
    if not matches:
        return pairs
    first, last = matches[0], matches[-1]
    out = []
    for k, p in enumerate(pairs):
        if k in aligned and (k < first or k > last):
            out.append((p[0], p[1], None, None))
            out.append((None, None, p[2], p[3]))
        else:
            out.append(p)
    return out


def gapless_segments(pairs):
    """Maximal runs of aligned residues with a constant offset, each cut
    back to its first and last identical residues:
    ``[(construct_beg, construct_end, ref_beg, ref_end), ...]``.

    A substitution inside a run (L33H) stays inside it; one at the edge of
    a gap (the P and G that replace Q141-K142 beside an inserted tag) is
    reported as a difference but does not start a mapped stretch.
    """
    runs = []
    current = None
    for cpos, cres, rpos, rres in pairs:
        if cpos is None or rpos is None:
            current = None
            continue
        offset = rpos - cpos
        if current and current[0] == offset:
            current[1].append((cpos, rpos, cres == rres))
            continue
        current = [offset, [(cpos, rpos, cres == rres)]]
        runs.append(current)
    segments = []
    for _offset, members in runs:
        matched = [m for m in members if m[2]]
        if not matched:
            continue
        (c0, r0, _), (c1, r1, _) = matched[0], matched[-1]
        segments.append((c0, c1, r0, r1))
    return segments


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
                # Construct residues the reference does not cover, at the
                # ends of the segment: the junction with whatever the user
                # cloned on, not a mapping to anything.
                out.append((cpos, cres, None, None, DIF_CLONING))
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
        pairs = trim_flanks(align_pairs(piece, self.ref_sequence))
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


def auth_numbering(structure):
    """``{chain: {label_seq_id: auth_seq_num}}`` from a structure whose
    residues carry ``label_seq`` (after ``assign_label_seq_id``); the author
    numbering is what ``pdbx_auth_*`` items record, and it need not equal
    the construct numbering (Ahir3252: author 11 is construct 1)."""
    out = {}
    for chain in structure[0]:
        mapping = out.setdefault(chain.name, {})
        for res in chain:
            if res.label_seq is not None:
                mapping.setdefault(res.label_seq, res.seqid.num)
    return out


def _auth(auth_map, chain, label, fallback=None):
    """The author number for a construct position; the nearest modelled
    neighbour's offset when the position itself is not modelled."""
    mapping = (auth_map or {}).get(chain) or {}
    if label in mapping:
        return mapping[label]
    if mapping:
        nearest = min(mapping, key=lambda k: abs(k - label))
        return mapping[nearest] + (label - nearest)
    return fallback if fallback is not None else label


def reference_tables(references, strands, entry_id="XXXX", extra_difs=(), auth_map=None):
    """The three categories as ``{category: (tags, rows)}``.

    ``references`` is a list of :class:`Reference` (any number per entity);
    ``strands`` maps entity id to its chain ids; ``extra_difs`` is
    ``[(entity_id, construct_pos, construct_res, details)]`` for residues no
    reference covers (tags, linkers), reported against the entity's first
    reference so OneDep sees them; ``auth_map`` (from :func:`auth_numbering`)
    gives the author numbering of each chain.
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
                                       ref.accession, r0, "?", r1, "?",
                                       _auth(auth_map, chain, c0), _auth(auth_map, chain, c1)))
            aid = first_align_for.get((ref.entity_id, chain))
            for cpos, cres, rpos, rres, details in ref.diffs:
                ordinal += 1
                struct_ref_seq_dif.append((
                    aid, entry_id, three_letter(cres) if cres else "?", chain,
                    cpos if cpos else "?", "?", ref.db_name, ref.accession,
                    three_letter(rres) if rres else "?", rpos if rpos else "?",
                    details, _auth(auth_map, chain, cpos) if cpos else "?", ordinal))
    for entity_id, cpos, cres, details in extra_difs:
        for chain in strands.get(str(entity_id), []):
            aid = first_align_for.get((str(entity_id), chain))
            if aid is None:
                continue
            ref = next(r for r in references if r.entity_id == str(entity_id))
            ordinal += 1
            struct_ref_seq_dif.append((
                aid, entry_id, three_letter(cres), chain, cpos, "?", ref.db_name,
                ref.accession, "?", "?", details, _auth(auth_map, chain, cpos), ordinal))
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
