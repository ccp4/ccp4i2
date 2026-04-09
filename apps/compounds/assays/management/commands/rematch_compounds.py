"""
Management command to retroactively match unmatched compound names in DataSeries.

Finds all DataSeries records where compound_name is set but compound FK is NULL,
and attempts to match them using the improved resolve_compound() utility.

Usage:
    # Dry run (preview matches without saving)
    python manage.py rematch_compounds --dry-run

    # Verbose dry run (show all unmatched names)
    python manage.py rematch_compounds --dry-run --verbose

    # Actually perform the matching
    python manage.py rematch_compounds

    # Match only for specific protocol(s)
    python manage.py rematch_compounds --protocol "ADME Microsome"

    # Match only for specific target(s)
    python manage.py rematch_compounds --target "Kinase X"
"""

import logging
from collections import defaultdict

from django.core.management.base import BaseCommand
from django.db import transaction

from compounds.assays.models import DataSeries
from compounds.utils import resolve_compound, resolve_compound_batch

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'Retroactively match unmatched compound names in DataSeries records'

    def add_arguments(self, parser):
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Preview matches without saving changes'
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show detailed output including unmatched names'
        )
        parser.add_argument(
            '--protocol',
            action='append',
            dest='protocols',
            help='Filter by protocol name (can be specified multiple times)'
        )
        parser.add_argument(
            '--target',
            action='append',
            dest='targets',
            help='Filter by target name (can be specified multiple times)'
        )
        parser.add_argument(
            '--batch-size',
            type=int,
            default=1000,
            help='Number of records to process per batch (default: 1000)'
        )

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        verbose = options['verbose']
        protocols = options.get('protocols') or []
        targets = options.get('targets') or []
        batch_size = options['batch_size']

        # Build queryset for unmatched DataSeries
        queryset = DataSeries.objects.filter(
            compound__isnull=True,
            compound_name__isnull=False,
        ).exclude(
            compound_name=''
        ).select_related(
            'assay',
            'assay__protocol',
            'assay__target',
        )

        # Apply filters
        if protocols:
            queryset = queryset.filter(assay__protocol__name__in=protocols)
        if targets:
            queryset = queryset.filter(
                assay__target__name__in=targets
            ) | queryset.filter(
                compound__target__name__in=targets
            )

        total_count = queryset.count()
        self.stdout.write(f"\nFound {total_count} DataSeries records with unmatched compound names\n")

        if total_count == 0:
            self.stdout.write(self.style.SUCCESS("No unmatched records to process."))
            return

        # Collect statistics
        matched_count = 0
        unmatched_count = 0
        matched_by_strategy = defaultdict(int)
        unmatched_names = defaultdict(int)
        updated_series = []

        # Process in batches for efficiency
        self.stdout.write(f"Processing in batches of {batch_size}...\n")

        # Collect all unique compound names first
        compound_names = list(queryset.values_list('compound_name', flat=True).distinct())
        self.stdout.write(f"Found {len(compound_names)} unique compound names to resolve\n")

        # Batch resolve all names
        resolution_map = resolve_compound_batch(compound_names)

        # Count matches
        for name, compound in resolution_map.items():
            if compound:
                matched_count += 1
            else:
                unmatched_count += 1
                unmatched_names[name] += 1

        self.stdout.write(f"\nResolution summary:")
        self.stdout.write(f"  - Matched: {matched_count} unique names")
        self.stdout.write(f"  - Unmatched: {unmatched_count} unique names\n")

        if dry_run:
            self.stdout.write(self.style.WARNING("\n=== DRY RUN - No changes will be saved ===\n"))

        # Now update DataSeries records
        records_updated = 0
        records_already_matched = 0

        with transaction.atomic():
            for series in queryset.iterator(chunk_size=batch_size):
                compound = resolution_map.get(series.compound_name)
                if compound:
                    if not dry_run:
                        series.compound = compound
                        series.save(update_fields=['compound'])
                    records_updated += 1
                    updated_series.append({
                        'series_id': str(series.id),
                        'compound_name': series.compound_name,
                        'matched_compound': compound.formatted_id,
                        'protocol': series.assay.protocol.name if series.assay and series.assay.protocol else 'N/A',
                    })

            if dry_run:
                # Rollback in dry run
                transaction.set_rollback(True)

        # Output results
        self.stdout.write(f"\n{'Would update' if dry_run else 'Updated'} {records_updated} DataSeries records\n")

        if verbose and updated_series:
            self.stdout.write("\nMatched records:")
            for item in updated_series[:50]:  # Limit output
                self.stdout.write(
                    f"  {item['compound_name']} -> {item['matched_compound']} "
                    f"(Protocol: {item['protocol']})"
                )
            if len(updated_series) > 50:
                self.stdout.write(f"  ... and {len(updated_series) - 50} more")

        if verbose and unmatched_names:
            self.stdout.write("\nStill unmatched compound names (top 50 by frequency):")
            sorted_unmatched = sorted(
                unmatched_names.items(),
                key=lambda x: -x[1]
            )[:50]

            # Import for diagnostic checks
            import re
            from compounds.registry.models import Compound
            # Note: malformed pattern must be tried FIRST to avoid incorrect matches
            ncl_malformed = re.compile(r'NCL0+-(\d+)', re.IGNORECASE)
            ncl_pattern = re.compile(r'NCL[-_\s]?0*(\d+)', re.IGNORECASE)

            for name, count in sorted_unmatched:
                # Count how many DataSeries have this name
                series_count = queryset.filter(compound_name=name).count()

                # Diagnose why NCL-format names didn't match
                # Try malformed pattern FIRST to avoid 0* consuming zeros incorrectly
                diagnosis = ""
                match = ncl_malformed.search(name) or ncl_pattern.search(name)
                if match:
                    reg_num = int(match.group(1))
                    exists = Compound.objects.filter(reg_number=reg_num).exists()
                    if not exists:
                        diagnosis = f" [reg_number {reg_num} NOT IN DATABASE]"
                    else:
                        diagnosis = f" [UNEXPECTED: reg_number {reg_num} exists!]"

                self.stdout.write(f"  '{name}' - {series_count} DataSeries records{diagnosis}")

        # Final summary
        self.stdout.write("\n" + "=" * 60)
        if dry_run:
            self.stdout.write(self.style.WARNING(
                f"DRY RUN COMPLETE: Would match {records_updated} of {total_count} records"
            ))
            self.stdout.write("Run without --dry-run to apply changes.")
        else:
            self.stdout.write(self.style.SUCCESS(
                f"COMPLETE: Matched {records_updated} of {total_count} records"
            ))

        remaining = total_count - records_updated
        if remaining > 0:
            self.stdout.write(f"\n{remaining} records remain unmatched.")
            self.stdout.write(
                "Consider adding aliases to compounds for frequently unmatched names."
            )
