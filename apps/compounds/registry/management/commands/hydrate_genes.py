"""
Hydrate Gene rows from HGNC.

For every Gene whose hydration_status is 'pending' or whose hydrated_at
is older than a configurable staleness window, call HGNC and fill in
name, aliases, uniprot_ids, ensembl_gene_id, hgnc_id.

Idempotent; safe to run on a schedule.

Usage:
    ccp4-python manage.py hydrate_genes
    ccp4-python manage.py hydrate_genes --dry-run
    ccp4-python manage.py hydrate_genes --symbol EGFR
    ccp4-python manage.py hydrate_genes --stale-days 30
    ccp4-python manage.py hydrate_genes --only-pending
"""

from datetime import timedelta

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone

from compounds.registry import hgnc
from compounds.registry.models import Gene


DEFAULT_STALE_DAYS = 90


class Command(BaseCommand):
    help = "Hydrate Gene rows from HGNC (aliases, cross-refs, canonical name)."

    def add_arguments(self, parser):
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Fetch from HGNC but do not write changes.',
        )
        parser.add_argument(
            '--symbol',
            type=str,
            help='Hydrate only the Gene with this symbol (case-insensitive).',
        )
        parser.add_argument(
            '--only-pending',
            action='store_true',
            help='Only hydrate Genes with status=pending; skip stale-but-ok ones.',
        )
        parser.add_argument(
            '--stale-days',
            type=int,
            default=DEFAULT_STALE_DAYS,
            help=f'Rehydrate Genes older than N days (default: {DEFAULT_STALE_DAYS}).',
        )
        parser.add_argument(
            '--timeout',
            type=int,
            default=hgnc.DEFAULT_TIMEOUT,
            help='HGNC request timeout in seconds.',
        )

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        only_symbol = options.get('symbol')
        only_pending = options['only_pending']
        stale_days = options['stale_days']
        timeout = options['timeout']

        queryset = self._build_queryset(only_symbol, only_pending, stale_days)
        total = queryset.count()

        if total == 0:
            self.stdout.write(self.style.SUCCESS('Nothing to hydrate.'))
            return

        self.stdout.write(f'Hydrating {total} Gene row(s)...')
        if dry_run:
            self.stdout.write(self.style.WARNING('DRY RUN — no changes will be saved.'))

        stats = {'ok': 0, 'not_found': 0, 'errors': 0}
        now = timezone.now()

        for gene in queryset.iterator():
            try:
                record = hgnc.fetch_by_symbol(gene.symbol, timeout=timeout)
            except hgnc.HGNCError as exc:
                stats['errors'] += 1
                self.stderr.write(
                    self.style.ERROR(f'  [{gene.symbol}] HGNC error: {exc}')
                )
                if not dry_run:
                    gene.hydration_status = 'failed'
                    gene.hydrated_at = now
                    gene.hydration_source = 'hgnc-api'
                    gene.save(update_fields=[
                        'hydration_status', 'hydrated_at', 'hydration_source'
                    ])
                continue

            if record is None:
                stats['not_found'] += 1
                self.stdout.write(
                    self.style.WARNING(f'  [{gene.symbol}] not found in HGNC')
                )
                if not dry_run:
                    gene.hydration_status = 'failed'
                    gene.hydrated_at = now
                    gene.hydration_source = 'hgnc-api'
                    gene.save(update_fields=[
                        'hydration_status', 'hydrated_at', 'hydration_source'
                    ])
                continue

            stats['ok'] += 1
            self._report_record(gene, record)

            if not dry_run:
                with transaction.atomic():
                    gene.symbol = record.symbol or gene.symbol  # canonicalise
                    gene.hgnc_id = record.hgnc_id
                    gene.name = record.name
                    gene.aliases = record.aliases
                    gene.uniprot_ids = record.uniprot_ids
                    gene.ensembl_gene_id = record.ensembl_gene_id
                    gene.hydration_status = 'ok'
                    gene.hydrated_at = now
                    gene.hydration_source = 'hgnc-api'
                    gene.save()

        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('Summary:'))
        self.stdout.write(f'  Hydrated:   {stats["ok"]}')
        self.stdout.write(f'  Not found:  {stats["not_found"]}')
        self.stdout.write(f'  Errors:     {stats["errors"]}')
        if dry_run:
            self.stdout.write(self.style.WARNING('DRY RUN — no changes were saved.'))

    # ---------------------------------------------------------------- #

    def _build_queryset(self, only_symbol, only_pending, stale_days):
        qs = Gene.objects.all()

        if only_symbol:
            qs = qs.filter(symbol__iexact=only_symbol.strip())
            if not qs.exists():
                raise CommandError(f'No Gene with symbol {only_symbol!r}')
            return qs

        if only_pending:
            return qs.filter(hydration_status='pending')

        # Default: pending OR (ok but stale)
        cutoff = timezone.now() - timedelta(days=stale_days)
        from django.db.models import Q
        return qs.filter(
            Q(hydration_status='pending') |
            Q(hydration_status='ok', hydrated_at__lt=cutoff) |
            Q(hydration_status='ok', hydrated_at__isnull=True) |
            Q(hydration_status='failed')  # retry previously-failed rows each run
        )

    def _report_record(self, gene, record):
        alias_preview = ', '.join(record.aliases[:3])
        extra = f' (+{len(record.aliases) - 3} more)' if len(record.aliases) > 3 else ''
        self.stdout.write(
            f'  [{gene.symbol}] {record.hgnc_id or "-"} '
            f'name="{record.name}" '
            f'aliases=[{alias_preview}{extra}]'
        )
