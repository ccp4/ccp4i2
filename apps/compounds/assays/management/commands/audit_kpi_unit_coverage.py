"""
Read-only audit of kpi_unit coverage across AnalysisResults.

Informs the NLP v1 scope decision (Q9 in apps/compounds/docs/NLP_QUERY_PROPOSAL.md):
can v1 drop the `Protocol.import_type='raw_data'` restriction, given that
`results['kpi_unit']` is the canonical unit source?

A row is "fully queryable" by NLP v1 when ALL of:
    - AnalysisResult.status == 'valid'
    - results['KPI'] is a non-empty string
    - results[results['KPI']] is numeric
    - results['kpi_unit'] is present
    - normalize_unit(results['kpi_unit']) is in VALID_UNITS

Usage:
    python manage.py audit_kpi_unit_coverage
    python manage.py audit_kpi_unit_coverage --verbose
    python manage.py audit_kpi_unit_coverage --json
    python manage.py audit_kpi_unit_coverage --protocol <uuid>

To run against a deployed Container App:
    az containerapp exec \\
        --name ccp4i2-bicep-server \\
        --resource-group "$RESOURCE_GROUP" \\
        --command "ccp4-python manage.py audit_kpi_unit_coverage"
"""

from collections import Counter, defaultdict

from django.core.management.base import BaseCommand

from compounds.assays.kpi_utils import (
    VALID_UNITS,
    normalize_unit,
    parse_unit_from_field_name,
)
from compounds.assays.models import AnalysisResult


class Command(BaseCommand):
    help = (
        "Audit kpi_unit coverage across AnalysisResults. Read-only. "
        "Informs the NLP v1 scope decision about import_type restriction."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Include sample problematic rows per bucket.',
        )
        parser.add_argument(
            '--json',
            action='store_true',
            help='Emit machine-readable JSON instead of the text report.',
        )
        parser.add_argument(
            '--protocol',
            type=str,
            help='Filter audit to a single Protocol UUID.',
        )
        parser.add_argument(
            '--chunk-size',
            type=int,
            default=1000,
            help='Queryset iteration chunk size (default: 1000).',
        )

    def handle(self, *args, **options):
        verbose = options['verbose']
        as_json = options['json']
        protocol_id = options.get('protocol')
        chunk_size = options['chunk_size']

        queryset = AnalysisResult.objects.select_related(
            'data_series__assay__protocol',
            'data_series__dilution_series',
        )
        if protocol_id:
            queryset = queryset.filter(data_series__assay__protocol_id=protocol_id)

        buckets = defaultdict(self._empty_bucket)

        for ar in queryset.iterator(chunk_size=chunk_size):
            ds = getattr(ar, 'data_series', None)
            protocol = ds.assay.protocol if (ds and ds.assay) else None
            import_type = protocol.import_type if protocol else '__no_protocol__'
            b = buckets[import_type]

            b['total'] += 1

            if ar.status != 'valid':
                continue
            b['valid_status'] += 1

            results = ar.results or {}
            kpi = results.get('KPI')
            if not isinstance(kpi, str) or not kpi:
                continue
            b['kpi_present'] += 1
            b['kpi_names'][kpi] += 1

            kpi_value = results.get(kpi)
            if not isinstance(kpi_value, (int, float)):
                continue
            b['kpi_numeric'] += 1

            kpi_unit = results.get('kpi_unit')
            if not kpi_unit:
                # Not currently queryable, but check whether populate_kpi_units
                # could recover the unit with its existing strategies.
                if ds and ds.dilution_series and ds.dilution_series.unit:
                    b['recoverable_from_ds'] += 1
                elif parse_unit_from_field_name(kpi):
                    b['recoverable_from_name'] += 1
                else:
                    b['unrecoverable'] += 1
                if verbose and len(b['sample_missing_unit']) < 5:
                    b['sample_missing_unit'].append(
                        f"id={ar.id} KPI={kpi!r} value={kpi_value}"
                    )
                continue
            b['kpi_unit_present'] += 1

            normalized = normalize_unit(kpi_unit) if isinstance(kpi_unit, str) else None
            if normalized not in VALID_UNITS:
                b['unknown_units'][str(kpi_unit)] += 1
                if verbose and len(b['sample_unknown_unit']) < 5:
                    b['sample_unknown_unit'].append(
                        f"id={ar.id} KPI={kpi!r} kpi_unit={kpi_unit!r}"
                    )
                continue
            b['kpi_unit_valid'] += 1
            b['fully_queryable'] += 1

        if as_json:
            self._emit_json(buckets)
        else:
            self._emit_text(buckets, verbose)

    # ------------------------------------------------------------------ #
    # helpers
    # ------------------------------------------------------------------ #

    def _empty_bucket(self):
        return {
            'total': 0,
            'valid_status': 0,
            'kpi_present': 0,
            'kpi_numeric': 0,
            'kpi_unit_present': 0,
            'kpi_unit_valid': 0,
            'fully_queryable': 0,
            'recoverable_from_ds': 0,
            'recoverable_from_name': 0,
            'unrecoverable': 0,
            'kpi_names': Counter(),
            'unknown_units': Counter(),
            'sample_missing_unit': [],
            'sample_unknown_unit': [],
        }

    def _emit_json(self, buckets):
        import json
        out = {}
        for imp, b in buckets.items():
            out[imp] = {
                **{k: v for k, v in b.items()
                   if k not in ('kpi_names', 'unknown_units',
                                'sample_missing_unit', 'sample_unknown_unit')},
                'kpi_names_top10': dict(b['kpi_names'].most_common(10)),
                'unknown_units_top10': dict(b['unknown_units'].most_common(10)),
                'sample_missing_unit': b['sample_missing_unit'],
                'sample_unknown_unit': b['sample_unknown_unit'],
            }
        self.stdout.write(json.dumps(out, indent=2, default=str))

    def _emit_text(self, buckets, verbose):
        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('kpi_unit coverage audit'))
        self.stdout.write('=' * 78)

        header = (
            f"{'import_type':<22}"
            f"{'total':>9}{'valid':>9}{'+KPI':>9}"
            f"{'+num':>9}{'+unit':>9}{'+valid':>9}{'ready%':>9}"
        )
        self.stdout.write(header)
        self.stdout.write('-' * len(header))

        totals = self._empty_bucket()
        totals['kpi_names'] = Counter()
        totals['unknown_units'] = Counter()

        for imp in sorted(buckets.keys()):
            b = buckets[imp]
            # ready% is the fraction of addressable rows (numeric-valued valid
            # rows) that are fully queryable by the NLP feature. This is the
            # relevant denominator for Q9: NLP filters out non-valid and
            # non-numeric rows anyway, so including them in the denominator
            # understates the real coverage.
            pct = (b['fully_queryable'] / b['kpi_numeric'] * 100.0) if b['kpi_numeric'] else 0.0
            row = (
                f"{imp:<22}"
                f"{b['total']:>9}{b['valid_status']:>9}{b['kpi_present']:>9}"
                f"{b['kpi_numeric']:>9}{b['kpi_unit_present']:>9}{b['kpi_unit_valid']:>9}"
                f"{pct:>8.1f}%"
            )
            if pct >= 80:
                self.stdout.write(self.style.SUCCESS(row))
            elif pct >= 50:
                self.stdout.write(self.style.WARNING(row))
            else:
                self.stdout.write(row)

            for key in ('total', 'valid_status', 'kpi_present', 'kpi_numeric',
                        'kpi_unit_present', 'kpi_unit_valid', 'fully_queryable',
                        'recoverable_from_ds', 'recoverable_from_name', 'unrecoverable'):
                totals[key] += b[key]
            totals['kpi_names'].update(b['kpi_names'])
            totals['unknown_units'].update(b['unknown_units'])

        # totals row
        self.stdout.write('-' * len(header))
        pct = (totals['fully_queryable'] / totals['kpi_numeric'] * 100.0) if totals['kpi_numeric'] else 0.0
        self.stdout.write(
            f"{'ALL':<22}"
            f"{totals['total']:>9}{totals['valid_status']:>9}{totals['kpi_present']:>9}"
            f"{totals['kpi_numeric']:>9}{totals['kpi_unit_present']:>9}{totals['kpi_unit_valid']:>9}"
            f"{pct:>8.1f}%"
        )

        self.stdout.write('')
        self.stdout.write('Column key:')
        self.stdout.write('  total   all AnalysisResult rows in bucket')
        self.stdout.write('  valid   + status == "valid"')
        self.stdout.write('  +KPI    + results["KPI"] is a non-empty string')
        self.stdout.write('  +num    + results[KPI] is numeric')
        self.stdout.write('  +unit   + results["kpi_unit"] is populated')
        self.stdout.write('  +valid  + normalize_unit(kpi_unit) is in VALID_UNITS')
        self.stdout.write('  ready%  fully_queryable / +num, % (NLP-addressable denominator)')

        # per-bucket detail
        for imp in sorted(buckets.keys()):
            b = buckets[imp]
            if b['total'] == 0:
                continue
            self.stdout.write('')
            self.stdout.write(self.style.SUCCESS(f'[{imp}]'))

            missing_unit_numeric = b['kpi_numeric'] - b['kpi_unit_present']
            if missing_unit_numeric > 0:
                self.stdout.write(
                    f'  {missing_unit_numeric} numeric-KPI rows are missing kpi_unit:'
                )
                self.stdout.write(
                    f'    recoverable via DilutionSeries: {b["recoverable_from_ds"]}'
                )
                self.stdout.write(
                    f'    recoverable via KPI field name: {b["recoverable_from_name"]}'
                )
                self.stdout.write(
                    f'    not recoverable: {b["unrecoverable"]}  '
                    f'(populate_kpi_units would not help these)'
                )

            if b['unknown_units']:
                self.stdout.write('  Unrecognised unit strings (top 10):')
                for unit, count in b['unknown_units'].most_common(10):
                    self.stdout.write(f'    {unit!r}: {count}')

            if b['kpi_names']:
                self.stdout.write('  KPI distribution (top 10):')
                for kpi, count in b['kpi_names'].most_common(10):
                    self.stdout.write(f'    {kpi!r}: {count}')

            if verbose:
                if b['sample_missing_unit']:
                    self.stdout.write('  Sample rows with missing kpi_unit:')
                    for s in b['sample_missing_unit']:
                        self.stdout.write(f'    {s}')
                if b['sample_unknown_unit']:
                    self.stdout.write('  Sample rows with unrecognised kpi_unit:')
                    for s in b['sample_unknown_unit']:
                        self.stdout.write(f'    {s}')

        # interpretation summary
        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('Interpretation for NLP v1 scope (Q9):'))
        ready_pct = (totals['fully_queryable'] / totals['kpi_numeric'] * 100.0) if totals['kpi_numeric'] else 0.0
        if ready_pct >= 90:
            self.stdout.write(
                f'  Overall addressable-row ready% = {ready_pct:.1f} — broad '
                'import_type coverage; scope restriction unnecessary.'
            )
        elif ready_pct >= 70:
            self.stdout.write(
                f'  Overall addressable-row ready% = {ready_pct:.1f} — most '
                'data is queryable; remaining gap should be surfaced via '
                'footer-count rather than scope restriction.'
            )
        else:
            self.stdout.write(
                f'  Overall addressable-row ready% = {ready_pct:.1f} — '
                'consider running populate_kpi_units and tightening import '
                'workflows before broadening NLP scope.'
            )
