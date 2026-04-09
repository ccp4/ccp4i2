"""
Management command to populate kpi_unit in AnalysisResult.results for existing data.

This command backfills unit information for existing analysis results by:
1. Copying unit from DilutionSeries (for dose-response fitted data)
2. Parsing unit from KPI field name (for table-of-values imports)

Usage:
    # Preview mode (dry run)
    python manage.py populate_kpi_units --dry-run

    # Run for all data
    python manage.py populate_kpi_units

    # Run for specific protocol
    python manage.py populate_kpi_units --protocol <protocol-id>

    # Show verbose output
    python manage.py populate_kpi_units --verbose
"""

from django.core.management.base import BaseCommand
from django.db import transaction

from compounds.assays.models import AnalysisResult, DataSeries
from compounds.assays.kpi_utils import normalize_unit, parse_unit_from_field_name


class Command(BaseCommand):
    """Populate kpi_unit in AnalysisResult.results for existing data."""

    help = "Backfill kpi_unit in analysis results from DilutionSeries or field names"

    def add_arguments(self, parser):
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Preview changes without saving',
        )
        parser.add_argument(
            '--protocol',
            type=str,
            help='Only process data series for a specific protocol UUID',
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show detailed output for each update',
        )
        parser.add_argument(
            '--batch-size',
            type=int,
            default=1000,
            help='Number of records to process in each batch (default: 1000)',
        )

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        protocol_id = options.get('protocol')
        verbose = options['verbose']
        batch_size = options['batch_size']

        if dry_run:
            self.stdout.write(self.style.WARNING('DRY RUN - no changes will be saved'))

        # Build queryset
        queryset = DataSeries.objects.select_related(
            'analysis',
            'dilution_series',
            'assay__protocol',
            'assay__protocol__preferred_dilutions',
        ).filter(
            analysis__isnull=False,
        )

        if protocol_id:
            queryset = queryset.filter(assay__protocol_id=protocol_id)
            self.stdout.write(f'Filtering to protocol: {protocol_id}')

        total_count = queryset.count()
        self.stdout.write(f'Processing {total_count} data series...')

        stats = {
            'already_has_unit': 0,
            'from_dilution_series': 0,
            'from_field_name': 0,
            'no_unit_found': 0,
            'errors': 0,
        }

        updated_results = []

        for i, ds in enumerate(queryset.iterator(chunk_size=batch_size)):
            if (i + 1) % 1000 == 0:
                self.stdout.write(f'  Processed {i + 1}/{total_count}...')

            try:
                analysis = ds.analysis
                if not analysis or not analysis.results:
                    stats['no_unit_found'] += 1
                    continue

                results = analysis.results

                # Skip if already has kpi_unit
                if results.get('kpi_unit'):
                    stats['already_has_unit'] += 1
                    continue

                unit = None
                source = None

                # Strategy 1: From DilutionSeries
                if ds.dilution_series and ds.dilution_series.unit:
                    unit = ds.dilution_series.unit
                    source = 'dilution_series'
                # Strategy 1b: From protocol's preferred_dilutions
                elif (ds.assay and ds.assay.protocol and
                      ds.assay.protocol.preferred_dilutions and
                      ds.assay.protocol.preferred_dilutions.unit):
                    unit = ds.assay.protocol.preferred_dilutions.unit
                    source = 'preferred_dilutions'
                # Strategy 2: Parse from KPI field name
                elif kpi_key := results.get('KPI'):
                    unit = parse_unit_from_field_name(kpi_key)
                    if unit:
                        source = 'field_name'

                if unit:
                    # Normalize the unit
                    unit = normalize_unit(unit)

                    if verbose:
                        self.stdout.write(
                            f'  {ds.id}: {results.get("KPI", "?")} -> {unit} (from {source})'
                        )

                    # Update results dict
                    results['kpi_unit'] = unit
                    analysis.results = results
                    updated_results.append(analysis)

                    if source in ('dilution_series', 'preferred_dilutions'):
                        stats['from_dilution_series'] += 1
                    else:
                        stats['from_field_name'] += 1
                else:
                    stats['no_unit_found'] += 1
                    if verbose:
                        self.stdout.write(
                            self.style.WARNING(
                                f'  {ds.id}: No unit found for KPI "{results.get("KPI", "?")}"'
                            )
                        )

            except Exception as e:
                stats['errors'] += 1
                self.stdout.write(
                    self.style.ERROR(f'  Error processing {ds.id}: {e}')
                )

        # Bulk update
        if updated_results and not dry_run:
            self.stdout.write(f'Saving {len(updated_results)} updates...')
            with transaction.atomic():
                AnalysisResult.objects.bulk_update(
                    updated_results,
                    ['results'],
                    batch_size=batch_size,
                )

        # Summary
        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('Summary:'))
        self.stdout.write(f'  Already had unit:      {stats["already_has_unit"]}')
        self.stdout.write(f'  From DilutionSeries:   {stats["from_dilution_series"]}')
        self.stdout.write(f'  From field name:       {stats["from_field_name"]}')
        self.stdout.write(f'  No unit found:         {stats["no_unit_found"]}')
        self.stdout.write(f'  Errors:                {stats["errors"]}')
        self.stdout.write(f'  Total updated:         {len(updated_results)}')

        if dry_run:
            self.stdout.write(self.style.WARNING('DRY RUN - no changes were saved'))
        else:
            self.stdout.write(self.style.SUCCESS('Done!'))
