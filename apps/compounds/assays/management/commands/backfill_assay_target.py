"""
Backfill Assay.target for assays that have target=NULL but whose DataSeries
compounds share an unambiguous target.

Assays created before the upload form auto-derived target (and any CLI imports)
have Assay.target=NULL. The aggregation page filters protocols by
`assays__target_id`, so those assays never surface.

Usage:
    # Preview what would change
    python manage.py backfill_assay_target --dry-run

    # Verbose preview (show per-assay reasoning)
    python manage.py backfill_assay_target --dry-run --verbose

    # Apply
    python manage.py backfill_assay_target

    # Limit to a protocol
    python manage.py backfill_assay_target --protocol "NCU LM"
"""

from collections import Counter

from django.core.management.base import BaseCommand
from django.db import transaction

from compounds.assays.models import Assay


class Command(BaseCommand):
    help = "Backfill Assay.target from matched compounds where unambiguous"

    def add_arguments(self, parser):
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Preview changes without saving',
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show per-assay reasoning',
        )
        parser.add_argument(
            '--protocol',
            action='append',
            dest='protocols',
            help='Limit to assays of this protocol (repeatable)',
        )

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        verbose = options['verbose']
        protocols = options.get('protocols') or []

        qs = Assay.objects.filter(target__isnull=True).select_related('protocol')
        if protocols:
            qs = qs.filter(protocol__name__in=protocols)

        total = qs.count()
        if total == 0:
            self.stdout.write("No assays with target=NULL found.")
            return

        self.stdout.write(f"Examining {total} assays with target=NULL"
                          + (" [DRY RUN]" if dry_run else ""))

        updated = 0
        ambiguous = 0
        no_compounds = 0
        ambiguous_details = []

        for assay in qs.iterator():
            target_ids = list(
                assay.data_series
                .filter(compound__isnull=False, compound__target__isnull=False)
                .values_list('compound__target_id', flat=True)
            )

            if not target_ids:
                no_compounds += 1
                if verbose:
                    self.stdout.write(
                        f"  SKIP  {assay.id} ({assay.protocol.name}): "
                        f"no matched compounds with target"
                    )
                continue

            counts = Counter(target_ids)
            if len(counts) > 1:
                ambiguous += 1
                ambiguous_details.append((assay, counts))
                if verbose:
                    breakdown = ", ".join(f"{tid}:{n}" for tid, n in counts.most_common())
                    self.stdout.write(
                        f"  SKIP  {assay.id} ({assay.protocol.name}): "
                        f"ambiguous targets [{breakdown}]"
                    )
                continue

            (sole_target_id,) = counts.keys()
            if verbose:
                self.stdout.write(
                    f"  SET   {assay.id} ({assay.protocol.name}) -> target {sole_target_id} "
                    f"(from {counts[sole_target_id]} series)"
                )

            if not dry_run:
                with transaction.atomic():
                    Assay.objects.filter(pk=assay.pk).update(target_id=sole_target_id)
            updated += 1

        self.stdout.write("")
        self.stdout.write("Summary:")
        self.stdout.write(f"  Updated:        {updated}"
                          + (" (would be)" if dry_run else ""))
        self.stdout.write(f"  Ambiguous:      {ambiguous}")
        self.stdout.write(f"  No compounds:   {no_compounds}")
        self.stdout.write(f"  Total examined: {total}")

        if ambiguous and not verbose:
            self.stdout.write("")
            self.stdout.write("Ambiguous assays (re-run with --verbose for full breakdown):")
            for assay, counts in ambiguous_details[:10]:
                breakdown = ", ".join(f"{tid}:{n}" for tid, n in counts.most_common())
                self.stdout.write(f"  {assay.id} ({assay.protocol.name}): [{breakdown}]")
            if len(ambiguous_details) > 10:
                self.stdout.write(f"  ...and {len(ambiguous_details) - 10} more")
