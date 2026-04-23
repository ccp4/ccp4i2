"""Protocol resolver tests (§6.2 + scope-limited pool).

The fixture shapes a realistic AR-programme protocol catalog plus an AKT
protocol that should NOT surface when scoped to AR, and a cross-project
selectivity assay (AR compounds tested against AKT) that exercises the
reg=X, assay=Y scope combination.
"""

from __future__ import annotations

import pytest

from compounds.assays.models import Assay, DataSeries, Protocol
from compounds.nlp.resolver import (
    resolve_protocol,
    resolve_targets,
    tokenize,
)
from compounds.nlp.spec import (
    FIELD_PROTOCOL_HINT,
    ProtocolClarify,
    ProtocolMiss,
    QuerySpec,
    ResolvedProtocol,
    ResolvedTargets,
)
from compounds.registry.models import Compound, Target


def _mk_ds(assay, compound, row=0):
    return DataSeries.objects.create(
        assay=assay,
        compound=compound,
        row=row,
        start_column=0,
        end_column=10,
    )


@pytest.fixture
def ar_activity(db):
    ar = Target.objects.create(name="AR degraders")
    akt = Target.objects.create(name="AKT programme")

    ar_cpd_1 = Compound.objects.create(target=ar, smiles="CCO")
    ar_cpd_2 = Compound.objects.create(target=ar, smiles="CCN")
    akt_cpd = Compound.objects.create(target=akt, smiles="CCC")

    # AR protocols (all on AR compounds, AR target)
    p_binding = Protocol.objects.create(name="AR binding HTRF")
    p_degradation = Protocol.objects.create(name="AR degradation HTRF")
    p_counter = Protocol.objects.create(name="AR counter-screen HTRF")
    p_v7 = Protocol.objects.create(name="AR-V7 IC50 TR-FRET")

    # AKT protocol (NOT on AR compounds) — should be filtered out
    p_akt = Protocol.objects.create(name="AKT inhibition HTRF")

    a_binding = Assay.objects.create(protocol=p_binding, target=ar)
    _mk_ds(a_binding, ar_cpd_1)
    _mk_ds(a_binding, ar_cpd_2, row=1)

    a_degradation = Assay.objects.create(protocol=p_degradation, target=ar)
    _mk_ds(a_degradation, ar_cpd_1)

    a_counter = Assay.objects.create(protocol=p_counter, target=ar)
    _mk_ds(a_counter, ar_cpd_1)

    a_v7 = Assay.objects.create(protocol=p_v7, target=ar)
    _mk_ds(a_v7, ar_cpd_1)

    a_akt = Assay.objects.create(protocol=p_akt, target=akt)
    _mk_ds(a_akt, akt_cpd)

    # Cross-project selectivity: AR compound run against AKT target
    p_selectivity = Protocol.objects.create(name="AKT selectivity HTRF")
    a_selectivity = Assay.objects.create(protocol=p_selectivity, target=akt)
    _mk_ds(a_selectivity, ar_cpd_1)

    return {
        "ar": ar,
        "akt": akt,
        "ar_cpd_1": ar_cpd_1,
        "p_binding": p_binding,
        "p_degradation": p_degradation,
        "p_counter": p_counter,
        "p_v7": p_v7,
        "p_akt": p_akt,
        "p_selectivity": p_selectivity,
    }


def _rt_both(target):
    """Build ResolvedTargets with registration=assay=target (default scope)."""
    return ResolvedTargets(registration=target, assay=target, scope_kind="both_same")


# ---------------------------------------------------------------------------
# tokenize()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("AR binding HTRF", {"ar", "binding", "htrf"}),
        ("AR-V7 IC50 TR-FRET", {"ar", "v7", "ic50", "tr", "fret"}),
        ("", set()),
        ("   ", set()),
        ("HTRF", {"htrf"}),
    ],
)
def test_tokenize(raw, expected):
    assert tokenize(raw) == expected


# ---------------------------------------------------------------------------
# Happy paths
# ---------------------------------------------------------------------------


def test_single_token_ambiguous_triggers_clarify(ar_activity):
    out = resolve_protocol("HTRF", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ProtocolClarify)
    names = {c.name for c in out.candidates}
    assert names == {
        "AR binding HTRF",
        "AR degradation HTRF",
        "AR counter-screen HTRF",
    }
    assert out.field == FIELD_PROTOCOL_HINT


def test_multi_token_narrowing_resolves(ar_activity):
    out = resolve_protocol("AR binding HTRF", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ResolvedProtocol)
    assert out.protocol.pk == ar_activity["p_binding"].pk


def test_partial_query_that_uniquely_matches_resolves(ar_activity):
    # "binding" appears only in "AR binding HTRF" among AR protocols.
    out = resolve_protocol("binding", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ResolvedProtocol)
    assert out.protocol.pk == ar_activity["p_binding"].pk


def test_punctuation_and_case_are_normalised(ar_activity):
    # "tr-fret" / "TR FRET" / "tr.fret" all tokenize the same way.
    for query in ("TR-FRET", "tr fret", "TR_FRET", "tr/fret"):
        out = resolve_protocol(query, _rt_both(ar_activity["ar"]))
        assert isinstance(out, ResolvedProtocol), query
        assert out.protocol.pk == ar_activity["p_v7"].pk


def test_candidates_carry_stats(ar_activity):
    out = resolve_protocol("HTRF", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ProtocolClarify)
    binding = next(c for c in out.candidates if c.name == "AR binding HTRF")
    assert binding.n_runs == 1
    assert binding.n_compounds == 2  # two compounds tested in that assay
    assert binding.last_run is not None


# ---------------------------------------------------------------------------
# Miss paths
# ---------------------------------------------------------------------------


def test_empty_hint_is_miss(ar_activity):
    out = resolve_protocol("", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ProtocolMiss)
    assert out.suggestions == []


def test_no_matching_protocol_returns_suggestions_by_overlap(ar_activity):
    # "binding elisa" — 1 token ("binding") overlaps with AR binding HTRF.
    out = resolve_protocol("binding elisa", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ProtocolMiss)
    assert 1 <= len(out.suggestions) <= 5
    # Best suggestion should be the one with most overlapping tokens.
    assert out.suggestions[0].name == "AR binding HTRF"


def test_completely_unknown_query_returns_empty_suggestions(ar_activity):
    # No tokens overlap with any in-scope protocol name.
    out = resolve_protocol(
        "obscureXYZ nonsenseQRS",
        _rt_both(ar_activity["ar"]),
    )
    assert isinstance(out, ProtocolMiss)
    assert out.suggestions == []


# ---------------------------------------------------------------------------
# Scope filtering
# ---------------------------------------------------------------------------


def test_scope_excludes_protocols_outside_registration_target(ar_activity):
    # AKT inhibition HTRF only tests AKT compounds — shouldn't appear when
    # we're scoped to AR compounds.
    out = resolve_protocol("HTRF", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ProtocolClarify)
    names = {c.name for c in out.candidates}
    assert "AKT inhibition HTRF" not in names


def test_cross_project_scope_finds_selectivity_protocol(ar_activity):
    # reg=AR, assay=AKT → the AKT selectivity run on an AR compound.
    rt = ResolvedTargets(
        registration=ar_activity["ar"],
        assay=ar_activity["akt"],
        scope_kind="cross",
    )
    out = resolve_protocol("selectivity HTRF", rt)
    assert isinstance(out, ResolvedProtocol)
    assert out.protocol.pk == ar_activity["p_selectivity"].pk


def test_cross_project_scope_excludes_akt_only_protocol(ar_activity):
    # AKT inhibition HTRF runs on AKT compounds (not AR compounds), so
    # even with assay=AKT, the registration=AR filter removes it.
    rt = ResolvedTargets(
        registration=ar_activity["ar"],
        assay=ar_activity["akt"],
        scope_kind="cross",
    )
    out = resolve_protocol("HTRF", rt)
    # Only one protocol matches in this scope: p_selectivity.
    assert isinstance(out, ResolvedProtocol)
    assert out.protocol.pk == ar_activity["p_selectivity"].pk


def test_assay_only_scope_finds_akt_native_protocol(ar_activity):
    # reg=None, assay=AKT → any assay targeting AKT regardless of compound
    # registration. AR compound's selectivity run AND AKT compound's
    # inhibition run both qualify.
    rt = ResolvedTargets(
        registration=None,
        assay=ar_activity["akt"],
        scope_kind="assay_only",
    )
    out = resolve_protocol("HTRF", rt)
    assert isinstance(out, ProtocolClarify)
    names = {c.name for c in out.candidates}
    assert names == {"AKT inhibition HTRF", "AKT selectivity HTRF"}


def test_empty_scope_returns_miss(ar_activity):
    # Target with no associated compounds/assays → no in-scope protocols.
    empty_target = Target.objects.create(name="Dormant programme")
    rt = ResolvedTargets(
        registration=empty_target,
        assay=empty_target,
        scope_kind="both_same",
    )
    out = resolve_protocol("HTRF", rt)
    assert isinstance(out, ProtocolMiss)
    assert out.suggestions == []


# ---------------------------------------------------------------------------
# Clarify ranking: fewer extra tokens first
# ---------------------------------------------------------------------------


def test_clarify_ordering_fewer_extra_tokens_first(ar_activity):
    # Query "AR HTRF" matches all three AR ... HTRF protocols equally well
    # on required tokens (each has both ar and htrf). Extras differ:
    #   AR binding HTRF       -> 1 extra (binding)
    #   AR degradation HTRF   -> 1 extra (degradation)
    #   AR counter-screen HTRF -> 2 extras (counter, screen)
    # So counter-screen must rank last.
    out = resolve_protocol("AR HTRF", _rt_both(ar_activity["ar"]))
    assert isinstance(out, ProtocolClarify)
    assert len(out.candidates) == 3
    assert out.candidates[-1].name == "AR counter-screen HTRF"


# ---------------------------------------------------------------------------
# Integration with the two-target orchestrator
# ---------------------------------------------------------------------------


def test_resolve_targets_then_resolve_protocol_roundtrip(ar_activity):
    spec = QuerySpec(
        registration_target_as_typed="AR degraders",
        assay_target_as_typed="AR degraders",
        protocol_hint="AR binding HTRF",
    )
    rt = resolve_targets(spec)
    assert isinstance(rt, ResolvedTargets)
    out = resolve_protocol(spec.protocol_hint, rt)
    assert isinstance(out, ResolvedProtocol)
    assert out.protocol.pk == ar_activity["p_binding"].pk
