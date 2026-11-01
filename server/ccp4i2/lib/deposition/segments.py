"""Proposing the source segments of a construct, so the depositor edits a
description rather than writing one.

The termini are where expression constructs carry what is not the protein:
a His tag, a protease site and the remnant its cut leaves. Those become a
``tag`` segment; the rest is one segment whose kind is guessed from the
sequence's name and description (a "designed" or "de novo" protein, a DARPin
scaffold), else ``natural`` with its UniProt accession left for the user.
Chimeras (a DARPin fused to an enzyme) still need splitting by hand: the
proposal says so rather than guess.

:func:`find_references` searches RCSB's sequence service for PDB entities
like a segment and returns their UniProt accessions with identity and
coverage, for the user to choose from; nothing is ever filled in unasked.
"""
import re

# N-terminal: up to six residues, a His run, then up to 25 residues to a
# protease site and the remnant its cut leaves (thrombin LVPR|GS, often
# followed by HM from an NdeI site; TEV ENLYFQ|G/S; HRV 3C LEVLFQ|GP).
_HIS_N = re.compile(r"^[A-Z]{0,6}?H{6,10}")
_PROTEASE = (
    ("thrombin", re.compile(r"LVPR[GS]S?(?:HM)?")),
    ("TEV", re.compile(r"ENLYFQ[GS]?(?:AM|HM)?")),
    ("HRV 3C", re.compile(r"LEVLFQ(?:GP)?")),
)
# C-terminal: an XhoI LE and/or a GS linker before six or more histidines.
_HIS_C = re.compile(r"(?:LE)?(?:GS)?H{6,10}$")

_DESIGNED_WORDS = re.compile(r"\b(designed|de ?novo|denovo|computationally designed|"
                             r"rfdiffusion|proteinmpnn|minibinder|binder design)\b", re.I)
_SCAFFOLD_WORDS = re.compile(r"\b(darpin|monobody|affibody|anticalin|adnectin|sso7d|"
                             r"nanobody|vhh|knottin|fynomer)\b", re.I)


def n_terminal_tag(sequence):
    """``(end, note)`` of an N-terminal tag (1-based, inclusive), or ``(0, None)``."""
    seq = sequence.upper()
    his = _HIS_N.match(seq)
    if not his:
        # A protease site right at the start (the tag was cut off but the
        # remnant of the site stayed): only a recognisable thrombin GSHM.
        if seq.startswith("GSHM"):
            return 4, "thrombin cleavage remnant GSHM"
        return 0, None
    end = his.end()
    window = seq[end:end + 25]
    for name, pattern in _PROTEASE:
        m = pattern.search(window)
        if m:
            return end + m.end(), "His tag and %s site, with the remnant of the cut" % name
    # No protease site: the His run and up to two G/S that close it.
    tail = re.match(r"[GS]{0,2}", seq[end:])
    return end + tail.end(), "His tag"


def c_terminal_tag(sequence):
    """``(start, note)`` of a C-terminal tag (1-based), or ``(None, None)``."""
    m = _HIS_C.search(sequence.upper())
    if not m:
        return None, None
    return m.start() + 1, "C-terminal His tag"


_GENERIC_WORDS = {"protein", "domain", "construct", "chain", "fusion", "the", "and", "with",
                  "binder", "complex", "variant", "mutant", "tag", "tagged"}


def _other_words(text, pattern):
    """Words of a name besides those ``pattern`` matches and generic ones:
    "Ahir-Darpin-3252" has "ahir" besides "darpin", so it is a fusion."""
    rest = pattern.sub(" ", text)
    return [w for w in re.findall(r"[A-Za-z]{3,}", rest) if w.lower() not in _GENERIC_WORDS]


def guess_kind(name="", description=""):
    """``designed``, ``scaffold``, ``natural``, or ``fusion`` (a scaffold or a
    design named together with another protein: proposed as natural, to be
    split by hand) from what the sequence is called."""
    text = " ".join(str(t or "") for t in (name, description))
    for pattern, kind in ((_SCAFFOLD_WORDS, "scaffold"), (_DESIGNED_WORDS, "designed")):
        if pattern.search(text):
            return "fusion" if _other_words(str(name or ""), pattern) else kind
    return "natural"


def propose_segments(sequence, name="", description=""):
    """``{"segments": [...], "notes": [...]}`` for one AU-contents sequence.

    Each segment is ``{"BEG", "END", "KIND"}`` (the def.xml field names, so
    the interface can write them as they are).
    """
    seq = "".join(str(sequence or "").split()).upper()
    n = len(seq)
    notes = []
    if n == 0:
        return {"segments": [], "notes": ["empty sequence"]}
    segments = []
    start = 1
    end_n, note_n = n_terminal_tag(seq)
    if end_n and end_n < n:
        segments.append({"BEG": 1, "END": end_n, "KIND": "tag"})
        notes.append("residues 1-%d look like %s" % (end_n, note_n))
        start = end_n + 1
    stop = n
    start_c, note_c = c_terminal_tag(seq)
    if start_c and start_c > start:
        stop = start_c - 1
    kind = guess_kind(name, description)
    if kind == "fusion":
        notes.append("the name suggests a fusion of a designed or scaffold protein with "
                     "another; split residues %d-%d into its parts" % (start, stop))
        kind = "natural"
    segments.append({"BEG": start, "END": stop, "KIND": kind})
    if kind == "natural":
        notes.append("residues %d-%d: give the UniProt accession (or use Find UniProt "
                     "reference), or change the kind if this is a design" % (start, stop))
    elif kind == "scaffold":
        notes.append("residues %d-%d: give the PDB entry of the scaffold" % (start, stop))
    if stop < n:
        segments.append({"BEG": stop + 1, "END": n, "KIND": "tag"})
        notes.append("residues %d-%d look like a %s" % (stop + 1, n, note_c))
    notes.append("a fusion of several proteins needs this segment split by hand")
    return {"segments": segments, "notes": notes}


RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/{entity}"


def find_references(sequence, max_hits=5, timeout=20, session=None):
    """PDB entities like ``sequence`` and their UniProt accessions, best
    first: ``{"ok": True, "hits": [{"pdb_entity", "uniprot", "identity",
    "coverage", "description"}]}``, or ``{"ok": False, "error"}`` when the
    service cannot be reached. Never raises."""
    try:
        import requests

        http = session or requests
        seq = "".join(str(sequence or "").split()).upper()
        if len(seq) < 20:
            return {"ok": False, "error": "sequence too short to search (fewer than 20 residues)"}
        query = {
            "query": {"type": "terminal", "service": "sequence",
                      "parameters": {"evalue_cutoff": 0.1, "identity_cutoff": 0.3,
                                     "sequence_type": "protein", "value": seq}},
            "return_type": "polymer_entity",
            "request_options": {"paginate": {"start": 0, "rows": max_hits * 3},
                                "results_verbosity": "verbose", "scoring_strategy": "sequence"},
        }
        response = http.post(RCSB_SEARCH, json=query, timeout=timeout)
        if response.status_code == 204:
            return {"ok": True, "hits": []}
        if response.status_code != 200:
            return {"ok": False, "error": "RCSB search: HTTP %d" % response.status_code}
        hits, seen = [], set()
        for result in response.json().get("result_set", []):
            identifier = result.get("identifier", "")
            identity = coverage = None
            for service in result.get("services", []):
                for node in service.get("nodes", []):
                    for ctx in node.get("match_context", []):
                        identity = ctx.get("sequence_identity", identity)
                        length = ctx.get("query_length") or len(seq)
                        aligned = ctx.get("query_aligned_seq")
                        if aligned:
                            coverage = len(aligned.replace("-", "")) / float(length)
            entry, _, entity = identifier.partition("_")
            detail = http.get(RCSB_ENTITY.format(entry=entry, entity=entity), timeout=timeout)
            if detail.status_code != 200:
                continue
            data = detail.json()
            ids = (data.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []
            description = (data.get("rcsb_polymer_entity") or {}).get("pdbx_description")
            for acc in ids:
                if acc in seen:
                    continue
                seen.add(acc)
                hits.append({"pdb_entity": identifier, "uniprot": acc, "identity": identity,
                             "coverage": round(coverage, 2) if coverage is not None else None,
                             "description": description})
            if len(hits) >= max_hits:
                break
        return {"ok": True, "hits": hits[:max_hits]}
    except Exception as err:  # offline, DNS, JSON: an answer, never a crash
        return {"ok": False, "error": "%s: %s" % (type(err).__name__, err)}


# Kinds whose segment needs a database reference to be described.
NEEDS_REFERENCE = {"natural": "UniProt accession", "engineered": "UniProt accession",
                   "scaffold": "PDB entry"}


def _runs(positions):
    runs = []
    for p in sorted(positions):
        if runs and runs[-1][1] == p - 1:
            runs[-1][1] = p
        else:
            runs.append([p, p])
    return [tuple(r) for r in runs]


def analyse_segments(sequence, segments):
    """How a list of segments covers a sequence, for display.

    ``segments`` are dicts with ``BEG``, ``END``, ``KIND`` and optionally
    ``UNIPROT``/``PDB_ID``. Returns ``{"length", "segments": [...],
    "uncovered": [(beg, end)], "overlaps": [(beg, end, [i, j])],
    "problems": [str]}``; each segment comes back with its index, length,
    the first and last residues, its reference ("none" when it has none)
    and its own problems. Never raises.
    """
    seq = "".join(str(sequence or "").split()).upper()
    n = len(seq)
    rows, problems = [], []
    owners = {}
    for i, seg in enumerate(segments or []):
        try:
            beg, end = int(seg.get("BEG") or 0), int(seg.get("END") or 0)
        except (TypeError, ValueError):
            beg, end = 0, 0
        kind = str(seg.get("KIND") or "natural")
        reference = (seg.get("UNIPROT") or seg.get("PDB_ID") or "").strip() or "none"
        own = []
        if beg < 1 or end < beg:
            own.append("range %s-%s is not a valid range" % (beg, end))
        if end > n:
            own.append("ends at %d, past the sequence end (%d)" % (end, n))
        if kind in NEEDS_REFERENCE and reference == "none":
            own.append("a %s segment needs a %s" % (kind, NEEDS_REFERENCE[kind]))
        lo, hi = max(beg, 1), min(end, n)
        for p in range(lo, hi + 1):
            owners.setdefault(p, []).append(i)
        piece = seq[lo - 1:hi] if hi >= lo else ""
        rows.append({"index": i, "BEG": beg, "END": end, "KIND": kind,
                     "length": max(0, end - beg + 1), "first": piece[:5], "last": piece[-5:],
                     "reference": reference, "problems": own})
        problems.extend("segment %d (%d-%d %s): %s" % (i + 1, beg, end, kind, p) for p in own)
    uncovered = _runs(p for p in range(1, n + 1) if p not in owners)
    overlapping = {}
    for p, who in owners.items():
        if len(who) > 1:
            overlapping.setdefault(tuple(who), []).append(p)
    overlaps = []
    for who, positions in overlapping.items():
        for beg, end in _runs(positions):
            overlaps.append((beg, end, list(who)))
    overlaps.sort()
    if uncovered:
        problems.append("residues covered by no segment: " +
                        ", ".join("%d-%d" % r if r[0] != r[1] else "%d" % r[0] for r in uncovered))
    for beg, end, who in overlaps:
        problems.append("residues %d-%d are in segments %s" % (
            beg, end, " and ".join(str(w + 1) for w in who)))
    return {"length": n, "segments": rows, "uncovered": uncovered, "overlaps": overlaps,
            "problems": problems}


def modelled_positions(sequence, model_path, polymer_type="PROTEIN"):
    """Positions of ``sequence`` (1-based) modelled in at least one chain of
    the model that aligns to it (identity 30% or more). Never raises."""
    try:
        import gemmi

        from .entities import residue_kind
        from ccp4i2.wrappers.modelASUCheck.script.modelASUCheck_report import coverage_from_alignment

        kind, ptype = residue_kind(polymer_type)
        full = gemmi.expand_one_letter_sequence("".join(str(sequence).split()).upper(), kind)
        st = gemmi.read_structure(str(model_path))
        st.setup_entities()
        covered = set()
        for chain in st[0]:
            poly = chain.get_polymer()
            if len(poly) == 0:
                continue
            result = gemmi.align_sequence_to_polymer(full, poly, ptype, gemmi.AlignmentScoring())
            if result.calculate_identity() < 30:
                continue
            _n, modelled, mismatched, _runs_ = coverage_from_alignment(
                result.add_gaps(gemmi.one_letter_code(full), 1),
                result.add_gaps(gemmi.one_letter_code(poly.extract_sequence()), 2),
                result.match_string)
            covered |= modelled | mismatched
        return sorted(covered)
    except Exception:
        return None
