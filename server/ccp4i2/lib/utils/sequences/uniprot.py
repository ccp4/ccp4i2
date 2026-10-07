"""UniProt from a name: candidates for "Human CDK2", and the sequence for one.

An agent asked to solve "CDK4/cyclin D" had a protein name and no sequence,
and nothing in CCP4i2 turned one into the other (the app's UniProt fetch
needs an accession and runs in the browser). This is the server's own.

``search(text, organism=None)`` reads what was typed, deterministically:

- an accession (``P24941``) or entry name (``CDK2_HUMAN``) is looked up directly;
- an organism is pulled out of the common shapes ("Human CDK2", "CDK2 from
  human", "CDK2 (human)", "CDK2, Homo sapiens", "CDK2 human") from a table of
  the organisms crystallographers meet, each with its taxonomy id;
- what is left is a gene symbol (one short token: exact gene name) or words
  (protein name), with UniProt's free text as the fallback.

It returns how it read the text with the ranked candidates (reviewed first,
then exact gene matches, then the organism asked for), and never picks one:
"cyclin D" is three genes. An ``organism`` argument wins over one in the text.
Only the text and the organism go to UniProt, never a sequence.

``fetch(accession, residue_range=None)`` gives the sequence and its
provenance, optionally cut to the crystallised construct ("175-432").
"""
import json
import re
import urllib.parse
import urllib.request

SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
FIELDS = "accession,id,protein_name,gene_primary,organism_name,organism_id,length,reviewed"
TIMEOUT = 30

# UniProt's documented accession pattern, and the entry-name form
ACCESSION = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$")
ENTRY_NAME = re.compile(r"^[A-Z0-9]{1,10}_[A-Z0-9]{1,5}$")
GENE_SYMBOL = re.compile(r"^[A-Za-z][A-Za-z0-9.\-]{0,15}$")

# The organisms crystallographers meet, by the names people type: taxonomy id,
# scientific name. Longest synonyms are tried first.
ORGANISMS = {
    9606: ("Homo sapiens", ["human", "homo sapiens", "h. sapiens", "h sapiens", "hs"]),
    10090: ("Mus musculus", ["mouse", "murine", "mus musculus", "m. musculus"]),
    10116: ("Rattus norvegicus", ["rat", "rattus norvegicus", "r. norvegicus"]),
    9913: ("Bos taurus", ["bovine", "cow", "bos taurus", "b. taurus"]),
    9823: ("Sus scrofa", ["pig", "porcine", "sus scrofa"]),
    9031: ("Gallus gallus", ["chicken", "gallus gallus"]),
    8355: ("Xenopus laevis", ["xenopus", "xenopus laevis", "frog"]),
    7955: ("Danio rerio", ["zebrafish", "danio rerio"]),
    7227: ("Drosophila melanogaster", ["drosophila", "fruit fly", "fly", "d. melanogaster",
                                       "drosophila melanogaster"]),
    6239: ("Caenorhabditis elegans", ["c. elegans", "c elegans", "worm", "caenorhabditis elegans"]),
    559292: ("Saccharomyces cerevisiae", ["yeast", "budding yeast", "s. cerevisiae",
                                          "saccharomyces cerevisiae"]),
    284812: ("Schizosaccharomyces pombe", ["fission yeast", "s. pombe", "schizosaccharomyces pombe"]),
    83333: ("Escherichia coli (K12)", ["e. coli", "e coli", "escherichia coli", "ecoli"]),
    224308: ("Bacillus subtilis", ["b. subtilis", "bacillus subtilis"]),
    83332: ("Mycobacterium tuberculosis", ["m. tuberculosis", "mtb", "tb",
                                           "mycobacterium tuberculosis"]),
    93061: ("Staphylococcus aureus", ["s. aureus", "staphylococcus aureus"]),
    208964: ("Pseudomonas aeruginosa", ["p. aeruginosa", "pseudomonas aeruginosa"]),
    300852: ("Thermus thermophilus", ["thermus", "t. thermophilus", "thermus thermophilus"]),
    3702: ("Arabidopsis thaliana", ["arabidopsis", "a. thaliana", "arabidopsis thaliana"]),
    2697049: ("SARS-CoV-2", ["sars-cov-2", "sars cov 2", "covid"]),
    11676: ("HIV-1", ["hiv-1", "hiv1", "hiv"]),
    36329: ("Plasmodium falciparum", ["p. falciparum", "plasmodium falciparum", "malaria"]),
}

_SYNONYMS = sorted(((syn, taxid) for taxid, (_, syns) in ORGANISMS.items() for syn in syns),
                   key=lambda s: -len(s[0]))


class UniProtError(RuntimeError):
    pass


def _organism_by_name(text):
    """(taxid, scientific name) for an organism as typed, or None."""
    key = " ".join(text.lower().replace("_", " ").split())
    if key.isdigit():
        taxid = int(key)
        return taxid, ORGANISMS.get(taxid, (key, []))[0]
    for synonym, taxid in _SYNONYMS:
        if key == synonym:
            return taxid, ORGANISMS[taxid][0]
    return None


def _pull_organism(text):
    """(rest of the text, (taxid, name) or None): an organism said in one of
    the common shapes, at the start or the end."""
    text = " ".join(text.split())
    for synonym, taxid in _SYNONYMS:
        s = re.escape(synonym)
        shapes = [
            rf"^(?:{s})[\s\-:]+(?P<rest>.+)$",                       # Human CDK2, human-CDK2
            rf"^(?P<rest>.+?)\s+(?:from|in|of)\s+(?:{s})$",          # CDK2 from human
            rf"^(?P<rest>.+?)\s*\(\s*(?:{s})\s*\)$",                 # CDK2 (human)
            rf"^(?P<rest>.+?)\s*,\s*(?:{s})$",                       # CDK2, Homo sapiens
            rf"^(?P<rest>.+?)\s+(?:{s})$",                           # CDK2 human
        ]
        for shape in shapes:
            m = re.match(shape, text, flags=re.IGNORECASE)
            if m and m.group("rest").strip():
                return m.group("rest").strip(), (taxid, ORGANISMS[taxid][0])
    return text, None


def read_query(text, organism=None):
    """How the text is read: {"kind": accession | entry_name | gene |
    protein_name, "term", "organism": {"taxid", "name"} or None}."""
    text = (text or "").strip()
    if not text:
        raise UniProtError("nothing to search for")
    upper = text.upper()
    if ACCESSION.match(upper):
        return {"kind": "accession", "term": upper, "organism": None}
    if ENTRY_NAME.match(upper) and "_" in upper:
        return {"kind": "entry_name", "term": upper, "organism": None}
    rest, found = _pull_organism(text)
    if organism not in (None, ""):
        found = _organism_by_name(str(organism)) or (None, str(organism))
    org = {"taxid": found[0], "name": found[1]} if found else None
    kind = "gene" if GENE_SYMBOL.match(rest) and " " not in rest else "protein_name"
    return {"kind": kind, "term": rest, "organism": org}


def _get(url, opener):
    request = urllib.request.Request(url, headers={"Accept": "application/json",
                                                   "User-Agent": "CCP4i2 (ccp4.ac.uk)"})
    try:
        with (opener or urllib.request.urlopen)(request, timeout=TIMEOUT) as response:
            return json.loads(response.read().decode("utf-8"))
    except Exception as err:  # noqa: BLE001 - offline, refused, or not JSON: say which
        raise UniProtError(f"UniProt could not be reached or answered oddly: {err}") from None


def _organism_clause(org):
    if not org:
        return ""
    if org.get("taxid"):
        return f" AND (organism_id:{org['taxid']})"
    return f' AND (organism_name:"{org["name"]}")'


def _query(reading, loose=False):
    term = reading["term"]
    if reading["kind"] == "gene" and not loose:
        core = f"(gene_exact:{term})"
    elif reading["kind"] == "protein_name" and not loose:
        core = f'(protein_name:"{term}")'
    else:
        core = f"({term})"
    return core + _organism_clause(reading.get("organism"))


def _candidate(entry):
    names = entry.get("proteinDescription", {})
    full = (names.get("recommendedName") or {}).get("fullName", {}).get("value")
    if not full:
        sub = names.get("submissionNames") or [{}]
        full = sub[0].get("fullName", {}).get("value")
    genes = entry.get("genes") or [{}]
    organism = entry.get("organism") or {}
    return {
        "accession": entry.get("primaryAccession"),
        "entry_name": entry.get("uniProtkbId"),
        "protein_name": full,
        "gene": (genes[0].get("geneName") or {}).get("value"),
        "organism": organism.get("scientificName"),
        "taxid": organism.get("taxonId"),
        "length": (entry.get("sequence") or {}).get("length"),
        "reviewed": "reviewed" in str(entry.get("entryType", "")).lower()
                    and "unreviewed" not in str(entry.get("entryType", "")).lower(),
    }


def _norm(text):
    """A protein name for comparing: lower case, hyphens and slashes as spaces."""
    return " ".join(re.sub(r"[-_/,]", " ", (text or "").lower()).split())


def _name_match(name, term):
    """How well a protein name matches the words typed, best 0: the same
    name; the name ending in the term and a short suffix ("G1/S-specific
    cyclin-D1" for "cyclin D"); the term as words in it; anything else."""
    name, term = _norm(name), _norm(term)
    if not term:
        return 3
    if name == term:
        return 0
    if re.search(rf"(^| ){re.escape(term)} ?[a-z0-9]{{0,2}}$", name):
        return 1
    if re.search(rf"(^| ){re.escape(term)}( |$)", name):
        return 2
    return 3


def _rank(candidates, reading):
    term = reading["term"]
    taxid = (reading.get("organism") or {}).get("taxid")

    def key(c):
        if reading["kind"] == "protein_name":
            match = _name_match(c["protein_name"], term)
        else:
            match = 0 if (c["gene"] or "").upper() == term.upper() else 3
        return (match > 1, not c["reviewed"], match,
                taxid is not None and c["taxid"] != taxid)
    return sorted(candidates, key=key)


def search(text, organism=None, limit=10, opener=None):
    """{"read_as": the reading, "query": what was asked of UniProt,
    "candidates": [...]} best first, at most ``limit``."""
    reading = read_query(text, organism)
    if reading["kind"] in ("accession", "entry_name"):
        entry = _get(ENTRY_URL.format(urllib.parse.quote(reading["term"])), opener)
        return {"read_as": reading, "query": reading["term"], "candidates": [_candidate(entry)]}
    queries = [_query(reading)]
    if reading["kind"] == "protein_name":
        # UniProt indexes "cyclin-D1" as one token, so the phrase "cyclin D"
        # misses the D-type cyclins: the words with the last as a prefix
        # finds them, and the ranking puts them first
        words = re.sub(r"[-/]", " ", reading["term"]).split()
        queries.append("protein_name:(" + " AND ".join(words[:-1] + [words[-1] + "*"]) + ")"
                       + _organism_clause(reading.get("organism")))
    queries.append(_query(reading, loose=True))
    found, asked = {}, []
    for i, query in enumerate(queries):
        if found and i == len(queries) - 1:
            break  # the loose free-text query only when the precise ones found nothing
        asked.append(query)
        url = SEARCH_URL + "?" + urllib.parse.urlencode(
            {"query": query, "fields": FIELDS, "format": "json", "size": max(limit * 5, 50)})
        for entry in _get(url, opener).get("results", []):
            candidate = _candidate(entry)
            found.setdefault(candidate["accession"], candidate)
    return {"read_as": reading, "query": " | ".join(asked),
            "candidates": _rank(list(found.values()), reading)[:limit]}


def parse_range(residue_range):
    """(first, last) from "175-432", 1-based inclusive; None for none."""
    if residue_range in (None, ""):
        return None
    m = re.match(r"^\s*(\d+)\s*[-:]\s*(\d+)\s*$", str(residue_range))
    if not m or int(m.group(1)) < 1 or int(m.group(2)) < int(m.group(1)):
        raise UniProtError(f"a residue range is 'first-last', e.g. 175-432, not {residue_range!r}")
    return int(m.group(1)), int(m.group(2))


def fetch(accession, residue_range=None, opener=None):
    """The entry's sequence and provenance, cut to ``residue_range`` if given:
    {"accession", "entry_name", "protein_name", "gene", "organism",
    "reviewed", "range", "sequence", "fasta"}."""
    key = (accession or "").strip().upper()
    if not (ACCESSION.match(key) or ENTRY_NAME.match(key)):
        raise UniProtError(f"{accession!r} is not a UniProt accession or entry name "
                           "(search by name first)")
    entry = _get(ENTRY_URL.format(urllib.parse.quote(key)), opener)
    sequence = (entry.get("sequence") or {}).get("value") or ""
    if not sequence:
        raise UniProtError(f"UniProt has no sequence for {key}")
    cut = parse_range(residue_range)
    if cut:
        first, last = cut
        if last > len(sequence):
            raise UniProtError(f"{key} has {len(sequence)} residues; the range ends at {last}")
        sequence = sequence[first - 1:last]
    out = _candidate(entry)
    out.update(range=f"{cut[0]}-{cut[1]}" if cut else None, sequence=sequence)
    db = "sp" if out["reviewed"] else "tr"
    span = f" residues {cut[0]}-{cut[1]}" if cut else ""
    out["fasta"] = (f">{db}|{out['accession']}|{out['entry_name']} {out['protein_name'] or ''}"
                    f" OS={out['organism'] or ''}{span}\n"
                    + "\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60)) + "\n")
    return out
