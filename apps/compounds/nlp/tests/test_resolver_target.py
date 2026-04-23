"""Target-field resolver tests.

The fixture mirrors the shape of DDU Database's hydrated target catalog
(post TARGET_MODEL_PROPOSAL.md backfill): a mix of single-gene targets,
multi-gene PPI targets, shared-gene pan-family ambiguity, and legacy
name-only fallbacks with no Gene links.
"""

from __future__ import annotations

import pytest

from compounds.nlp.resolver import normalize, resolve_target
from compounds.nlp.spec import ResolvedTarget, TargetClarify, TargetMiss
from compounds.registry.models import Gene, Target


@pytest.fixture
def ddu_catalog(db):
    """Build a mini target catalog shaped like the DDU production data."""
    egfr = Gene.objects.create(
        symbol="EGFR",
        name="epidermal growth factor receptor",
        aliases=["ERBB", "ERBB1", "HER1"],
        hydration_status="ok",
    )
    skp2 = Gene.objects.create(symbol="SKP2", aliases=["FBL1", "FBXL1"], hydration_status="ok")
    cks1b = Gene.objects.create(symbol="CKS1B", aliases=["CKS1"], hydration_status="ok")
    aurka = Gene.objects.create(symbol="AURKA", aliases=["AIK", "ARK1"], hydration_status="ok")
    aurkb = Gene.objects.create(symbol="AURKB", aliases=["AIK2", "ARK2"], hydration_status="ok")
    myc = Gene.objects.create(symbol="MYC", aliases=["BHLHE39", "C-MYC"], hydration_status="ok")

    egfr_target = Target.objects.create(name="EGFR programme")
    egfr_target.genes.add(egfr)

    skp2_cks1 = Target.objects.create(name="Skp2-Cks1")
    skp2_cks1.genes.add(skp2, cks1b)

    myc_aur = Target.objects.create(name="Myc-Aur")
    myc_aur.genes.add(aurka, aurkb, myc)

    myc_rna = Target.objects.create(name="Myc RNA")
    myc_rna.genes.add(myc)

    ar_deg = Target.objects.create(name="AR degraders")  # legacy name-only

    return {
        "egfr": egfr_target,
        "skp2_cks1": skp2_cks1,
        "myc_aur": myc_aur,
        "myc_rna": myc_rna,
        "ar_degraders": ar_deg,
    }


# ---------------------------------------------------------------------------
# normalize() — pure function, no DB needed
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "raw,expected",
    [
        ("EGFR", "egfr"),
        ("  egfr  ", "egfr"),
        ("Skp2-Cks1", "skp2cks1"),
        ("AR degraders", "ardegraders"),
        ("ERBB1", "erbb1"),
        ("", ""),
        ("   ", ""),
    ],
)
def test_normalize(raw, expected):
    assert normalize(raw) == expected


# ---------------------------------------------------------------------------
# Exact matches across every pool-entry flavour
# ---------------------------------------------------------------------------

def test_resolve_exact_name_hit(ddu_catalog):
    res = resolve_target("EGFR programme")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["egfr"].pk
    assert res.matched_via == "name"


def test_resolve_exact_gene_symbol(ddu_catalog):
    res = resolve_target("EGFR")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["egfr"].pk
    assert res.matched_via == "gene_symbol"


def test_resolve_exact_gene_alias(ddu_catalog):
    res = resolve_target("ERBB1")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["egfr"].pk
    assert res.matched_via == "gene_alias"


def test_resolve_exact_gene_name(ddu_catalog):
    res = resolve_target("epidermal growth factor receptor")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["egfr"].pk
    assert res.matched_via == "gene_name"


# ---------------------------------------------------------------------------
# Normalisation — case, whitespace, punctuation
# ---------------------------------------------------------------------------

def test_resolve_case_insensitive(ddu_catalog):
    res = resolve_target("egfr")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["egfr"].pk


def test_resolve_whitespace_tolerant(ddu_catalog):
    res = resolve_target("  SKP2  ")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["skp2_cks1"].pk
    assert res.matched_via == "gene_symbol"


def test_resolve_hyphen_stripped_name(ddu_catalog):
    # "Skp2-Cks1" and "Skp2Cks1" both normalise to "skp2cks1" → hits the
    # Target.name pool entry for Skp2-Cks1.
    for query in ("Skp2-Cks1", "skp2 cks1", "Skp2Cks1"):
        res = resolve_target(query)
        assert isinstance(res, ResolvedTarget), query
        assert res.target.pk == ddu_catalog["skp2_cks1"].pk
        assert res.matched_via == "name"


def test_resolve_ppi_via_either_partner(ddu_catalog):
    for query in ("SKP2", "CKS1B", "CKS1"):  # the CKS1 alias
        res = resolve_target(query)
        assert isinstance(res, ResolvedTarget), query
        assert res.target.pk == ddu_catalog["skp2_cks1"].pk


# ---------------------------------------------------------------------------
# Name-only fallback (legacy, unhydrated)
# ---------------------------------------------------------------------------

def test_resolve_name_only_target(ddu_catalog):
    res = resolve_target("AR degraders")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["ar_degraders"].pk
    assert res.matched_via == "name"


def test_resolve_name_only_ignores_punctuation(ddu_catalog):
    res = resolve_target("ar_degraders")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["ar_degraders"].pk


# ---------------------------------------------------------------------------
# Ambiguity — shared gene across multiple targets
# ---------------------------------------------------------------------------

def test_resolve_shared_gene_triggers_clarify(ddu_catalog):
    res = resolve_target("MYC")
    assert isinstance(res, TargetClarify)
    ids = {c.id for c in res.candidates}
    assert ids == {str(ddu_catalog["myc_aur"].id), str(ddu_catalog["myc_rna"].id)}


def test_clarify_candidates_carry_gene_symbols(ddu_catalog):
    res = resolve_target("MYC")
    assert isinstance(res, TargetClarify)
    by_id = {c.id: c for c in res.candidates}
    myc_aur = by_id[str(ddu_catalog["myc_aur"].id)]
    assert myc_aur.gene_symbols == ["AURKA", "AURKB", "MYC"]


# ---------------------------------------------------------------------------
# Fuzzy match — distance-1, gated on ≥4-char normalised query
# ---------------------------------------------------------------------------

def test_transposition_is_not_distance_one(ddu_catalog):
    # "egrf" vs "egfr": same length, two positions differ (indices 1 and 2).
    # That is Levenshtein distance 2; transpositions are not treated as
    # single edits. So the typo is not auto-corrected — user gets a miss.
    res = resolve_target("EGRF")
    assert isinstance(res, TargetMiss)


def test_resolve_fuzzy_single_insertion(ddu_catalog):
    # "egfrs" vs "egfr" — one insertion, distance 1, query is 5 chars.
    res = resolve_target("egfrs")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["egfr"].pk
    assert res.matched_via == "gene_symbol_fuzzy"


def test_resolve_fuzzy_single_substitution(ddu_catalog):
    # "egfx" vs "egfr" — single substitution, distance 1, query 4 chars.
    res = resolve_target("egfx")
    assert isinstance(res, ResolvedTarget)
    assert res.target.pk == ddu_catalog["egfr"].pk
    assert res.matched_via == "gene_symbol_fuzzy"


def test_fuzzy_suppressed_below_threshold(ddu_catalog):
    # "myk" is distance 1 from gene symbol MYC but only 3 chars long —
    # §6.1 says fuzzy requires ≥4 chars on the normalised query.
    res = resolve_target("MYK")
    assert isinstance(res, TargetMiss)


def test_fuzzy_suppressed_for_short_ambiguous_codenames():
    # Without any fixtures this would fall through anyway; keep the
    # explicit assertion so §6.1's "prevents AR matching AKT/AXL/ABL"
    # rationale has a regression test even in isolation.
    assert len(normalize("AR")) < 4


# ---------------------------------------------------------------------------
# Miss — returns top-5 suggestions
# ---------------------------------------------------------------------------

def test_resolve_miss_returns_suggestions(ddu_catalog):
    res = resolve_target("notarealtarget")
    assert isinstance(res, TargetMiss)
    assert res.query == "notarealtarget"
    assert 1 <= len(res.suggestions) <= 5
    assert all(s.name for s in res.suggestions)


def test_resolve_empty_query_is_miss(ddu_catalog):
    res = resolve_target("")
    assert isinstance(res, TargetMiss)
    assert res.suggestions == []


def test_resolve_whitespace_only_query_is_miss(ddu_catalog):
    res = resolve_target("   ")
    assert isinstance(res, TargetMiss)
    assert res.suggestions == []


# ---------------------------------------------------------------------------
# Pool construction — an un-genned target still matches by name
# ---------------------------------------------------------------------------

def test_name_only_target_excluded_from_gene_queries(ddu_catalog):
    # "AR degraders" has no Gene links; querying the symbol "AR" should
    # not resolve to it (AR isn't in the pool at all, and normalised "ar"
    # is under the fuzzy threshold regardless).
    res = resolve_target("AR")
    assert isinstance(res, TargetMiss)
