"""Tanimoto similarity narrowing for the NLP "compounds similar to X"
feature (slice 21).

Anchor compounds are resolved via the existing ``resolve_compound_ref``
(slice 14). Each anchor's pre-computed Morgan fingerprint (radius=2,
2048 bits, stored on MolecularProperties via the registration-time
signal hook) is loaded once. The candidate compound set is filtered
in Python: a compound is kept if its fingerprint is ≥ the threshold
similar to ANY anchor (UNION semantics — the chemist asking *"similar
to A and B"* gets neighbours of either; intersection is harder to
phrase usefully and easier to add later).

Compounds without a fingerprint (registered before slice 21, or with
SMILES that didn't parse) are silently dropped from the result set —
the chemist gets the analogue of the *"compounds with no signal"*
case in slice 19's ranking. Run the ``compute_fingerprints`` management
command to backfill.
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Set


# De-facto medicinal-chemistry cutoff. Below 0.5 results get noisy;
# 0.7 is "looks like the same series"; ≥0.85 is essentially "near-
# duplicate or close analogue". Default chosen for "similar enough
# to be interesting to a SAR-thinking chemist".
DEFAULT_SIMILAR_THRESHOLD = 0.7


def _load_morgan_fp(bit_string: Optional[str]):
    """Deserialise a stored Morgan fingerprint (2048-char "0101..." bit
    string) back to an RDKit ExplicitBitVect, or None when the value
    is missing or unreadable. We use bit-string round-trip rather than
    ToBinary/CreateFromBinaryText because the latter doesn't preserve
    the bit-vector length on the round-trip — see model docstring."""
    if not bit_string:
        return None
    try:
        from rdkit import DataStructs
    except ImportError:
        return None
    try:
        return DataStructs.CreateFromBitString(bit_string)
    except Exception:
        return None


def neighbours_within_threshold(
    anchor_compounds: Iterable,
    candidate_compound_pks: Set,
    threshold: float,
) -> Set:
    """Return the subset of ``candidate_compound_pks`` whose Morgan
    fingerprint has Tanimoto similarity ≥ ``threshold`` to ANY of the
    anchor compounds. Anchors are always included in the result if
    they're in the candidate set.

    Defensive against missing fingerprints — a compound without one
    is silently excluded, same as a compound that doesn't pass the
    threshold. Use ``compute_fingerprints --recompute`` to backfill
    missing ones.
    """
    try:
        from rdkit import DataStructs
    except ImportError:
        return set()

    anchor_fps = []
    for c in anchor_compounds:
        props = getattr(c, "molecular_properties", None)
        if props is None:
            continue
        fp = _load_morgan_fp(props.morgan_fp)
        if fp is not None:
            anchor_fps.append(fp)
    if not anchor_fps:
        # No anchor has a usable fingerprint — nothing to compare to.
        # Don't fall through to "match everything"; that would be an
        # invisible miss. Return empty set.
        return set()

    # Lazy import — keep the registry module out of the import graph
    # for callers that never invoke similarity.
    from compounds.registry.models import MolecularProperties

    matching: Set = set()
    rows = (
        MolecularProperties.objects
        .filter(compound_id__in=candidate_compound_pks)
        .exclude(morgan_fp__isnull=True)
        .values_list("compound_id", "morgan_fp")
    )
    for compound_pk, fp_bytes in rows:
        candidate_fp = _load_morgan_fp(fp_bytes)
        if candidate_fp is None:
            continue
        # UNION across anchors — keep on first match. Bulk Tanimoto is
        # one C call per anchor; for typical anchor counts (1-3) a
        # short Python loop is fine.
        best = max(
            DataStructs.TanimotoSimilarity(candidate_fp, anchor_fp)
            for anchor_fp in anchor_fps
        )
        if best >= threshold:
            matching.add(compound_pk)
    return matching
