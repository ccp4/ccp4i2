"""
HGNC (HUGO Gene Nomenclature Committee) REST client.

Fetches gene metadata from https://rest.genenames.org. Used to hydrate
the `Gene` model with aliases, cross-references, and canonical names.

Pure stdlib; no new dependency. HGNC is free, unauthenticated, and
industry-standard for human gene nomenclature. See
apps/compounds/docs/TARGET_MODEL_PROPOSAL.md §5.1.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)

HGNC_REST_BASE = "https://rest.genenames.org"
USER_AGENT = "ccp4i2-compounds/1.0 (+https://github.com/ccp4/ccp4i2)"
DEFAULT_TIMEOUT = 10  # seconds


class HGNCError(Exception):
    """Raised when the HGNC API returns a non-200 or malformed response."""


@dataclass
class GeneRecord:
    """Normalized HGNC gene record."""

    symbol: str
    hgnc_id: str = ""
    name: str = ""
    aliases: list[str] = field(default_factory=list)
    uniprot_ids: list[str] = field(default_factory=list)
    ensembl_gene_id: str = ""


def fetch_by_symbol(symbol: str, timeout: int = DEFAULT_TIMEOUT) -> Optional[GeneRecord]:
    """
    Fetch HGNC metadata for the given approved gene symbol.

    Returns None if the symbol is not known to HGNC (empty docs list).
    Raises HGNCError on network/transport errors.
    """
    if not symbol or not symbol.strip():
        return None
    normalized = symbol.strip().upper()
    url = f"{HGNC_REST_BASE}/fetch/symbol/{quote(normalized)}"
    payload = _fetch_json(url, timeout)
    docs = _extract_docs(payload)
    if not docs:
        return None
    return _parse_doc(docs[0])


def search(query: str, timeout: int = DEFAULT_TIMEOUT) -> list[GeneRecord]:
    """
    Text-search HGNC for candidates matching the query (e.g. a target name).

    Returns a list of GeneRecord, empty if no match. The HGNC search ranks
    by relevance; we preserve that order.
    """
    if not query or not query.strip():
        return []
    url = f"{HGNC_REST_BASE}/search/{quote(query.strip())}"
    payload = _fetch_json(url, timeout)
    docs = _extract_docs(payload)
    # Search returns lightweight records; for full metadata the caller should
    # follow up with fetch_by_symbol on each hit's symbol. We still parse
    # what's present so callers can use the list directly when only symbol
    # + hgnc_id is enough.
    return [_parse_doc(d) for d in docs]


def _fetch_json(url: str, timeout: int) -> dict:
    req = Request(
        url,
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    try:
        with urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode("utf-8")
    except HTTPError as exc:
        raise HGNCError(f"HGNC returned HTTP {exc.code} for {url}") from exc
    except URLError as exc:
        raise HGNCError(f"HGNC network error for {url}: {exc.reason}") from exc
    try:
        return json.loads(body)
    except json.JSONDecodeError as exc:
        raise HGNCError(f"HGNC returned non-JSON body for {url}: {exc}") from exc


def _extract_docs(payload: dict) -> list[dict]:
    """HGNC wraps results as {responseHeader: {...}, response: {numFound, docs: [...]}}."""
    try:
        return list(payload["response"]["docs"])
    except (KeyError, TypeError):
        return []


def _parse_doc(doc: dict) -> GeneRecord:
    # HGNC keeps alias_symbol and prev_symbol as separate lists; we union them,
    # dedupe while preserving order, and drop empties.
    raw_aliases: list[str] = []
    for key in ("alias_symbol", "prev_symbol"):
        value = doc.get(key) or []
        if isinstance(value, list):
            raw_aliases.extend(v for v in value if isinstance(v, str))

    seen: set[str] = set()
    aliases: list[str] = []
    for alias in raw_aliases:
        alias_stripped = alias.strip()
        if alias_stripped and alias_stripped not in seen:
            seen.add(alias_stripped)
            aliases.append(alias_stripped)

    uniprot_ids_raw = doc.get("uniprot_ids") or []
    uniprot_ids = [u for u in uniprot_ids_raw if isinstance(u, str) and u.strip()]

    return GeneRecord(
        symbol=(doc.get("symbol") or "").strip(),
        hgnc_id=(doc.get("hgnc_id") or "").strip(),
        name=(doc.get("name") or "").strip(),
        aliases=aliases,
        uniprot_ids=uniprot_ids,
        ensembl_gene_id=(doc.get("ensembl_gene_id") or "").strip(),
    )
