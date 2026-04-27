"""Backfill Morgan fingerprints on MolecularProperties (slice 21).

The fingerprint is computed at registration time for new compounds (via
the post-save signal that creates MolecularProperties). This command
fills in fingerprints for compounds registered BEFORE the slice-21
migration, plus any rows where the fingerprint is null because RDKit
choked on a particular SMILES at the time.

Run via cron / Azure scheduled job, or manually after deploying the
migration:
    ccp4-python manage.py compute_fingerprints
    ccp4-python manage.py compute_fingerprints --dry-run
    ccp4-python manage.py compute_fingerprints --recompute   # replaces existing fp
"""

from __future__ import annotations

from django.core.management.base import BaseCommand
from django.db.models import Q

from compounds.registry.models import Compound, MolecularProperties


class Command(BaseCommand):
    help = "Compute Morgan fingerprints for MolecularProperties without one."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="List how many compounds need fingerprints without computing.",
        )
        parser.add_argument(
            "--recompute",
            action="store_true",
            help="Recompute fingerprints even when present (e.g. after an RDKit upgrade).",
        )
        parser.add_argument(
            "--limit",
            type=int,
            default=None,
            help="Process at most N compounds (useful for incremental backfills).",
        )

    def handle(
        self, *args,
        dry_run: bool = False, recompute: bool = False,
        limit: int | None = None, **options,
    ):
        if recompute:
            qs = Compound.objects.all()
        else:
            # Compounds with no MolecularProperties row OR with a row
            # where morgan_fp is still null.
            qs = Compound.objects.filter(
                Q(molecular_properties__isnull=True)
                | Q(molecular_properties__morgan_fp__isnull=True)
            ).distinct()

        total = qs.count()
        if limit is not None:
            qs = qs[:limit]

        n_to_process = qs.count()
        if n_to_process == 0:
            self.stdout.write("Nothing to do — every compound has a fingerprint.")
            return

        self.stdout.write(
            f"{n_to_process} compound(s) to process "
            f"(of {total} total{' missing fingerprint' if not recompute else ''})."
        )
        if dry_run:
            self.stdout.write("DRY RUN — no changes made.")
            return

        n_done = 0
        n_failed = 0
        for compound in qs.iterator():
            result = MolecularProperties.calculate_for_compound(compound)
            if result is None or result.morgan_fp is None:
                n_failed += 1
            else:
                n_done += 1
            if (n_done + n_failed) % 200 == 0:
                self.stdout.write(
                    f"  ... {n_done} computed, {n_failed} failed"
                )

        self.stdout.write(
            self.style.SUCCESS(
                f"Done: {n_done} fingerprint(s) computed, "
                f"{n_failed} failed (likely SMILES that don't parse)."
            )
        )
