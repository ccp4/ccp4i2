"""
Heuristically populate Protocol.target_value / poor_value / threshold_scale
for existing protocols that have no thresholds configured.

Three categories drive the defaults:

1. Plate-based potency (raw_data, ms_intact, or name matches SPR/Kd/Ki/ITC)
     target = 10 nM, poor = 10,000 nM (10 uM), log scale
     -> IC50/EC50/Kd-style binding assays with fitted values

2. Cellular (table_of_values, not a binding assay)
     target = 500 nM, poor = 30,000 nM (30 uM), log scale
     -> growth inhibition, HTRF, viability (pre-fit value tables)

3. ADME (pharmaron_adme or name matches an ADME keyword)
     per-endpoint defaults; see ADME_RULES below. Unmatched ADME
     protocols are left alone — hand-curate from the Protocol edit UI.

Defaults assume the KPI is reported in the expected units for each
assay type (nM for potency, uM for CYP inhibition, uL/min/mg for
microsomal clearance, etc.). Anomalous units on a given protocol will
produce nonsensical thresholds — review the output of --dry-run.

Usage:
    # Preview (nothing is written)
    python manage.py backfill_protocol_thresholds --dry-run

    # Apply to all unconfigured protocols
    python manage.py backfill_protocol_thresholds

    # Re-apply to everything including already-configured protocols
    python manage.py backfill_protocol_thresholds --overwrite

    # Limit to protocols whose name contains a fragment (repeatable)
    python manage.py backfill_protocol_thresholds --protocol Microsom --dry-run
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass

from django.core.management.base import BaseCommand
from django.db import transaction

from compounds.assays.models import Protocol


@dataclass
class Classification:
    target_value: float
    poor_value: float
    threshold_scale: str  # 'log' | 'linear'
    reason: str


# ADME sub-type keyword -> thresholds. First match wins; order matters
# because "solubility" is more specific than many generic words.
# Each pattern is compiled case-insensitive; use \b for word boundaries
# where a keyword could appear as a substring of a larger word.
ADME_RULES: list[tuple[re.Pattern[str], Classification | None]] = [
    # ========================================================================
    # Explicit skips (no sensible two-anchor default) — matched first so they
    # take precedence over the classification rules below.
    # ========================================================================

    # Physicochemical — wrong axis (pH) or band-around-midpoint shape.
    (re.compile(r'log\s*[dp](?![a-z])', re.IGNORECASE), None),
    (re.compile(r'pk[ab](?![a-z])', re.IGNORECASE), None),

    # In-vivo PK parameters: units vary wildly (mL/min/kg, L/kg, h, %),
    # cannot default safely. Matches on "in vivo", bioavailability,
    # half-life, AUC, Vd, "Clearance", or the unit strings mL/min/L/kg.
    (
        re.compile(
            r'\bin\s*vivo\b'
            r'|bioavail'
            r'|half[\s_-]?life'
            r'|(?<![a-z])auc(?![a-z])'
            r'|(?<![a-z])vd(?![a-z])'
            r'|clearance'
            r'|(?<![a-z])m[lL]\s*min'
            r'|\bl\s*kg\b',
            re.IGNORECASE,
        ),
        None,
    ),

    # Intact MS — KPI conventions not standardised across protocols; skip.
    (re.compile(r'(?<![a-z])intact(?![a-z])', re.IGNORECASE), None),

    # ========================================================================
    # Classification rules.
    # ========================================================================

    (
        re.compile(r'solubility|aq[_\s]*sol(?![a-z])|(?<![a-z])pbs(?![a-z])|(?<![a-z])sol(?![a-z])', re.IGNORECASE),
        Classification(100, 10, 'log', 'ADME: solubility (higher-better)'),
    ),
    (
        re.compile(r'microsom|(?<![a-z])[hmr]lm(?![a-z])', re.IGNORECASE),
        Classification(10, 50, 'log', 'ADME: microsomal clearance (lower-better)'),
    ),
    (
        re.compile(r'hepatocyt|(?<![a-z])[hmr]hep(?![a-z])', re.IGNORECASE),
        Classification(10, 100, 'log', 'ADME: hepatocyte clearance (lower-better)'),
    ),
    (
        re.compile(r'permeab|caco|mdck|papp', re.IGNORECASE),
        Classification(10, 1, 'log', 'ADME: permeability (higher-better)'),
    ),
    # Plasma protein binding: accept underscore or whitespace between
    # "plasma" and "protein" so Pharmaron's `plasma_protein_binding` matches.
    (
        re.compile(r'(?<![a-z])ppb(?![a-z])|plasma[_\s]+protein|%\s*free', re.IGNORECASE),
        Classification(20, 1, 'log', 'ADME: plasma protein binding (% free)'),
    ),
    # Matrix stability: GSH, plasma, blood, serum — all % compound remaining
    # after incubation. Use negative lookaround so "gsh_stability" matches
    # even though "_" is a word character (so \b fails).
    (
        re.compile(r'(?<![a-z])gsh(?![a-z])|glutathione|(?:blood|serum|plasma)[_\s]*stab', re.IGNORECASE),
        Classification(90, 50, 'linear', 'ADME: matrix stability (% remaining)'),
    ),
    (
        re.compile(r'(?<![a-z])cyp\d*(?![a-z])|cytochrome', re.IGNORECASE),
        Classification(30, 1, 'log', 'ADME: CYP inhibition (higher IC50 = better)'),
    ),
]

# Binding-assay name override — treat SPR/Kd/Ki/ITC as plate-based potency
# even when the data was imported as table_of_values.
# Uses negative lookaround (rather than \b) so `ITCinput` and similar
# non-word-boundary forms are still recognised.
BINDING_NAME = re.compile(r'(?<![a-z])(?:spr|k[di]|itc)(?![a-z])', re.IGNORECASE)

PLATE_POTENCY = Classification(10, 10_000, 'log', 'Plate-based potency')
CELLULAR = Classification(500, 30_000, 'log', 'Cellular (table_of_values)')


def classify(protocol: Protocol) -> Classification | None:
    """
    Return the threshold classification for a protocol, or None if the
    command should leave it alone (unclassifiable or explicitly skipped).
    """
    name = protocol.name or ''

    # ADME first — both because import_type flags it and because name
    # keywords can apply to non-ADME-typed protocols that happen to
    # measure an ADME-style endpoint.
    for pattern, classification in ADME_RULES:
        if pattern.search(name):
            return classification  # may be None (skip, e.g. logD)

    # pharmaron_adme with no recognised sub-type — bail, hand-curate
    if protocol.import_type == 'pharmaron_adme':
        return None

    # Intact MS — KPI conventions not standardised; leave unscored.
    if protocol.import_type == 'ms_intact':
        return None

    # SPR/Kd/Ki/ITC override: potency despite table_of_values import
    if BINDING_NAME.search(name):
        return Classification(
            PLATE_POTENCY.target_value,
            PLATE_POTENCY.poor_value,
            PLATE_POTENCY.threshold_scale,
            'Binding assay (SPR/Kd/Ki/ITC) — potency defaults',
        )

    if protocol.import_type == 'raw_data':
        return PLATE_POTENCY

    if protocol.import_type == 'table_of_values':
        return CELLULAR

    return None


class Command(BaseCommand):
    help = (
        "Heuristically populate Protocol.target_value / poor_value / "
        "threshold_scale for protocols missing thresholds"
    )

    def add_arguments(self, parser):
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Print what would change without writing',
        )
        parser.add_argument(
            '--overwrite',
            action='store_true',
            help='Also re-classify protocols that already have thresholds set',
        )
        parser.add_argument(
            '--protocol',
            action='append',
            dest='protocol_filters',
            default=[],
            help='Limit to protocols whose name contains this substring '
                 '(repeatable, case-insensitive)',
        )

    def handle(self, *args, **options):
        dry_run: bool = options['dry_run']
        overwrite: bool = options['overwrite']
        filters: list[str] = options['protocol_filters']

        qs = Protocol.objects.all().order_by('name')
        if filters:
            from django.db.models import Q
            q = Q()
            for f in filters:
                q |= Q(name__icontains=f)
            qs = qs.filter(q)

        if dry_run:
            self.stdout.write(self.style.WARNING(
                '\n=== DRY RUN — no changes will be saved ===\n'
            ))

        counts: Counter[str] = Counter()
        applied = 0
        preserved = 0
        unclassified: list[str] = []

        for protocol in qs:
            already_set = (
                protocol.target_value is not None
                and protocol.poor_value is not None
            )
            if already_set and not overwrite:
                counts['preserved (already set)'] += 1
                preserved += 1
                continue

            classification = classify(protocol)
            if classification is None:
                counts['unclassified'] += 1
                unclassified.append(
                    f'{protocol.name} [import_type={protocol.import_type}]'
                )
                continue

            counts[classification.reason] += 1

            action = 'WOULD SET' if dry_run else 'SET'
            overwrite_note = ' (overwriting)' if already_set else ''
            self.stdout.write(
                f'  {action}{overwrite_note}  {protocol.name}'
            )
            self.stdout.write(
                f'        target={classification.target_value}  '
                f'poor={classification.poor_value}  '
                f'scale={classification.threshold_scale}  '
                f'[{classification.reason}]'
            )

            if not dry_run:
                protocol.target_value = classification.target_value
                protocol.poor_value = classification.poor_value
                protocol.threshold_scale = classification.threshold_scale
                # Using transaction-per-save rather than bulk_update so that
                # any model-level clean() runs — there are few protocols and
                # this is a rare one-off, so performance is not a concern.
                with transaction.atomic():
                    protocol.save(update_fields=[
                        'target_value', 'poor_value', 'threshold_scale',
                    ])

            applied += 1

        # Summary
        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('=== Summary ==='))
        self.stdout.write(f'  Protocols scanned:  {qs.count()}')
        self.stdout.write(
            f'  {"Would apply" if dry_run else "Applied"}: {applied}'
        )
        self.stdout.write(f'  Preserved (already configured): {preserved}')
        self.stdout.write(f'  Unclassified (no change): {len(unclassified)}')
        self.stdout.write('')
        for reason, count in sorted(counts.items(), key=lambda x: -x[1]):
            self.stdout.write(f'    {count:4d}  {reason}')

        if unclassified:
            self.stdout.write('')
            self.stdout.write(self.style.WARNING(
                'Unclassified protocols (hand-curate from the Protocol edit UI):'
            ))
            for name in unclassified:
                self.stdout.write(f'    - {name}')

        if dry_run:
            self.stdout.write('')
            self.stdout.write(self.style.NOTICE(
                'Re-run without --dry-run to apply these changes.'
            ))
