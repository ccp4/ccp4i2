"""Morgan fingerprint integration tests (slice 21).

The fingerprint is computed alongside the descriptors (MW, cLogP, etc.)
in ``MolecularProperties.calculate_for_compound``, which fires from the
Compound.save() hook on every SMILES change. These tests pin that
behaviour and the backfill management command.
"""

from __future__ import annotations

import pytest
from django.contrib.auth import get_user_model
from django.core.management import call_command
from io import StringIO

from compounds.registry.models import Compound, MolecularProperties, Target

User = get_user_model()


# ---------------------------------------------------------------------------
# Fingerprint computed at registration
# ---------------------------------------------------------------------------


def test_fingerprint_computed_on_compound_save(db):
    target = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=target, smiles="c1ccncc1")  # pyridine
    props = MolecularProperties.objects.get(compound=c)
    assert props.morgan_fp is not None
    # Stored as a 2048-char bit string.
    assert len(props.morgan_fp) == 2048
    assert set(props.morgan_fp) <= {"0", "1"}


def test_fingerprint_recomputed_on_smiles_change(db):
    target = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=target, smiles="c1ccncc1")
    fp_pyridine = c.molecular_properties.morgan_fp

    c.smiles = "c1ccc2[nH]ncc2c1"  # indazole — completely different molecule
    c.save()
    c.molecular_properties.refresh_from_db()
    fp_indazole = c.molecular_properties.morgan_fp
    assert fp_pyridine != fp_indazole


def test_invalid_smiles_does_not_crash_save(db):
    """SMILES that RDKit can't parse should leave morgan_fp null
    rather than raising or losing the descriptor row."""
    target = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=target, smiles="not-a-valid-smiles")
    props = MolecularProperties.objects.filter(compound=c).first()
    if props is not None:
        # The whole props row may be skipped (calculate_for_compound
        # returns None on parse failure); if it does get created, the
        # fingerprint should be null.
        assert props.morgan_fp is None


# ---------------------------------------------------------------------------
# Tanimoto similarity helper
# ---------------------------------------------------------------------------


def test_tanimoto_self_similarity_is_one(db):
    """A compound is exactly similar to itself."""
    from compounds.nlp.similarity import neighbours_within_threshold
    target = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=target, smiles="c1ccncc1")
    matches = neighbours_within_threshold(
        anchor_compounds=[c],
        candidate_compound_pks={c.pk},
        threshold=0.99,
    )
    assert c.pk in matches


def test_unrelated_compounds_below_threshold(db):
    """Aspirin and caffeine have very low Tanimoto (~0.09)."""
    from compounds.nlp.similarity import neighbours_within_threshold
    target = Target.objects.create(name="ARd")
    aspirin = Compound.objects.create(target=target, smiles="CC(=O)Oc1ccccc1C(=O)O")
    caffeine = Compound.objects.create(target=target, smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    matches = neighbours_within_threshold(
        anchor_compounds=[aspirin],
        candidate_compound_pks={aspirin.pk, caffeine.pk},
        threshold=0.5,
    )
    assert aspirin.pk in matches
    assert caffeine.pk not in matches


def test_close_analogues_above_threshold(db):
    """Aspirin and salicylic acid are close analogues (~0.45 Tanimoto)."""
    from compounds.nlp.similarity import neighbours_within_threshold
    target = Target.objects.create(name="ARd")
    aspirin = Compound.objects.create(target=target, smiles="CC(=O)Oc1ccccc1C(=O)O")
    salicylic = Compound.objects.create(target=target, smiles="Oc1ccccc1C(=O)O")
    matches = neighbours_within_threshold(
        anchor_compounds=[aspirin],
        candidate_compound_pks={aspirin.pk, salicylic.pk},
        threshold=0.4,
    )
    assert salicylic.pk in matches


def test_anchor_with_no_fingerprint_returns_empty(db):
    """If the anchor's fingerprint is missing (e.g. backfill not run
    on a pre-existing compound), the helper returns empty rather
    than silently matching everything."""
    from compounds.nlp.similarity import neighbours_within_threshold
    target = Target.objects.create(name="ARd")
    pyridine = Compound.objects.create(target=target, smiles="c1ccncc1")
    # Strip the anchor's fingerprint to simulate the un-backfilled case.
    MolecularProperties.objects.filter(compound=pyridine).update(morgan_fp=None)
    pyridine.refresh_from_db()
    matches = neighbours_within_threshold(
        anchor_compounds=[pyridine],
        candidate_compound_pks={pyridine.pk},
        threshold=0.7,
    )
    assert matches == set()


def test_multi_anchor_union_keeps_neighbours_of_either(db):
    """UNION across anchors: a compound matching either anchor is kept.
    Use aspirin and caffeine as two unrelated anchor scaffolds; salicylic
    acid is close to aspirin (~0.45) and theophylline is close to
    caffeine (~0.46); both should land at threshold 0.4."""
    from compounds.nlp.similarity import neighbours_within_threshold
    target = Target.objects.create(name="ARd")
    aspirin = Compound.objects.create(target=target, smiles="CC(=O)Oc1ccccc1C(=O)O")
    caffeine = Compound.objects.create(
        target=target, smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    )
    near_aspirin = Compound.objects.create(target=target, smiles="Oc1ccccc1C(=O)O")
    near_caffeine = Compound.objects.create(
        target=target, smiles="Cn1c(=O)c2[nH]cnc2n(C)c1=O",
    )
    matches = neighbours_within_threshold(
        anchor_compounds=[aspirin, caffeine],
        candidate_compound_pks={
            aspirin.pk, caffeine.pk, near_aspirin.pk, near_caffeine.pk,
        },
        threshold=0.4,
    )
    assert near_aspirin.pk in matches
    assert near_caffeine.pk in matches


# ---------------------------------------------------------------------------
# Backfill management command
# ---------------------------------------------------------------------------


def test_backfill_command_dry_run_does_not_compute(db):
    target = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=target, smiles="c1ccncc1")
    MolecularProperties.objects.filter(compound=c).update(morgan_fp=None)
    out = StringIO()
    call_command("compute_fingerprints", "--dry-run", stdout=out)
    c.molecular_properties.refresh_from_db()
    assert c.molecular_properties.morgan_fp is None
    assert "DRY RUN" in out.getvalue()


def test_backfill_command_computes_missing_fingerprints(db):
    target = Target.objects.create(name="ARd")
    c1 = Compound.objects.create(target=target, smiles="c1ccncc1")
    c2 = Compound.objects.create(target=target, smiles="c1cncnc1")
    MolecularProperties.objects.filter(
        compound__in=[c1, c2],
    ).update(morgan_fp=None)
    out = StringIO()
    call_command("compute_fingerprints", stdout=out)
    c1.molecular_properties.refresh_from_db()
    c2.molecular_properties.refresh_from_db()
    assert c1.molecular_properties.morgan_fp is not None
    assert c2.molecular_properties.morgan_fp is not None


def test_backfill_command_idempotent_when_all_present(db):
    target = Target.objects.create(name="ARd")
    Compound.objects.create(target=target, smiles="c1ccncc1")
    out = StringIO()
    call_command("compute_fingerprints", stdout=out)
    assert "Nothing to do" in out.getvalue()


def test_backfill_command_recompute_replaces_existing(db):
    target = Target.objects.create(name="ARd")
    c = Compound.objects.create(target=target, smiles="c1ccncc1")
    original = c.molecular_properties.morgan_fp
    # Stub a different value to confirm --recompute overwrites it.
    MolecularProperties.objects.filter(compound=c).update(
        morgan_fp="0" * 2048,
    )
    call_command("compute_fingerprints", "--recompute", stdout=StringIO())
    c.molecular_properties.refresh_from_db()
    assert c.molecular_properties.morgan_fp == original
