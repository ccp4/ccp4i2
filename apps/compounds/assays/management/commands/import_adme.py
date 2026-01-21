"""
Management command to import ADME assay data from vendor files.

Supports auto-detection of file type and vendor based on filename patterns.
Currently supports NCU (National Crystallography Unit) ADME data:
- LM (Liver Microsome Stability)
- BS (Blood/Serum Stability)
- GSH (Glutathione Stability)
- Caco-2 (Permeability)

Usage:
    # Preview mode (dry run)
    python manage.py import_adme /path/to/ADME-NCU-LM-20231218.xlsx --preview

    # Import with compound matching
    python manage.py import_adme /path/to/ADME-NCU-LM-20231218.xlsx

    # Import multiple files
    python manage.py import_adme file1.xlsx file2.xlsx file3.xlsx

    # Force specific parser
    python manage.py import_adme /path/to/file.xlsx --parser ncu-lm
"""

import re
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from compounds.assays.importers import (
    detect_parser,
    get_parser,
    list_parsers,
    ParseResult,
)
from compounds.assays.models import Assay, DataSeries, AnalysisResult, Protocol
from compounds.registry.models import Compound


class Command(BaseCommand):
    """Import ADME assay data from vendor Excel files."""

    help = "Import ADME assay data from vendor files (NCU LM, BS, GSH, Caco-2)"

    def add_arguments(self, parser):
        parser.add_argument(
            'files',
            nargs='+',
            type=str,
            help='Path(s) to ADME data file(s)',
        )
        parser.add_argument(
            '--preview',
            action='store_true',
            help='Preview import without saving (dry run)',
        )
        parser.add_argument(
            '--parser',
            type=str,
            help='Force specific parser by protocol slug (e.g., ncu-lm)',
        )
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show detailed output',
        )
        parser.add_argument(
            '--skip-unmatched',
            action='store_true',
            help='Skip compounds that cannot be matched to registry',
        )
        parser.add_argument(
            '--list-parsers',
            action='store_true',
            help='List available parsers and exit',
        )

    def handle(self, *args, **options):
        # List parsers and exit if requested
        if options['list_parsers']:
            self._list_parsers()
            return

        files = options['files']
        preview = options['preview']
        parser_slug = options.get('parser')
        verbose = options['verbose']
        skip_unmatched = options['skip_unmatched']

        self.verbose = verbose

        if preview:
            self.stdout.write(self.style.WARNING(
                "\n[PREVIEW MODE - no changes will be saved]\n"
            ))

        total_results = 0
        total_errors = 0
        total_matched = 0
        total_unmatched = 0

        for filepath in files:
            path = Path(filepath)
            if not path.exists():
                self.stderr.write(self.style.ERROR(f"File not found: {filepath}"))
                total_errors += 1
                continue

            self.stdout.write(f"\n{'=' * 60}")
            self.stdout.write(f"Processing: {path.name}")
            self.stdout.write('=' * 60)

            try:
                result = self._process_file(
                    path,
                    parser_slug=parser_slug,
                    preview=preview,
                    skip_unmatched=skip_unmatched,
                )

                total_results += result['results_count']
                total_errors += result['errors_count']
                total_matched += result['matched_count']
                total_unmatched += result['unmatched_count']

            except Exception as e:
                self.stderr.write(self.style.ERROR(f"Failed to process {path.name}: {e}"))
                if verbose:
                    import traceback
                    self.stderr.write(traceback.format_exc())
                total_errors += 1

        # Summary
        self.stdout.write(f"\n{'=' * 60}")
        self.stdout.write("IMPORT SUMMARY")
        self.stdout.write('=' * 60)
        self.stdout.write(f"Files processed: {len(files)}")
        self.stdout.write(f"Total results: {total_results}")
        self.stdout.write(f"  Matched to compounds: {total_matched}")
        self.stdout.write(f"  Unmatched: {total_unmatched}")
        self.stdout.write(f"Errors: {total_errors}")

        if preview:
            self.stdout.write(self.style.WARNING(
                "\n[PREVIEW MODE - no changes were saved]"
            ))

    def _list_parsers(self):
        """List all available parsers."""
        parsers = list_parsers()
        self.stdout.write("\nAvailable ADME Parsers:")
        self.stdout.write("-" * 60)

        for slug, parser in sorted(parsers.items(), key=lambda x: x[0]):
            self.stdout.write(f"  {slug:<20} {parser.vendor} {parser.assay_type}")

        self.stdout.write(f"\nTotal: {len(parsers)} parsers")

    def _process_file(
        self,
        filepath: Path,
        parser_slug: str = None,
        preview: bool = False,
        skip_unmatched: bool = False,
    ) -> dict:
        """Process a single ADME file."""

        # Get parser
        if parser_slug:
            parser = get_parser(parser_slug)
            if not parser:
                raise CommandError(f"Unknown parser: {parser_slug}")
        else:
            parser = detect_parser(filepath)
            if not parser:
                raise CommandError(
                    f"Could not detect parser for: {filepath.name}\n"
                    "Use --parser to specify explicitly or --list-parsers to see options."
                )

        self.stdout.write(f"Using parser: {parser.protocol_slug} ({parser.vendor} {parser.assay_type})")

        # Parse file
        parse_result: ParseResult = parser.parse(filepath)

        # Report errors
        if parse_result.errors:
            self.stdout.write(self.style.WARNING(f"\nParsing warnings/errors:"))
            for error in parse_result.errors:
                style = self.style.ERROR if error.severity == 'error' else self.style.WARNING
                self.stdout.write(style(f"  [{error.severity}] {error.message}"))

        if not parse_result.results:
            self.stdout.write(self.style.WARNING("No results parsed from file"))
            return {
                'results_count': 0,
                'errors_count': len([e for e in parse_result.errors if e.severity == 'error']),
                'matched_count': 0,
                'unmatched_count': 0,
            }

        self.stdout.write(f"\nParsed {len(parse_result.results)} results")

        # Match compounds
        matched_count = 0
        unmatched_count = 0
        results_to_import = []

        for parsed_result in parse_result.results:
            compound_id = parsed_result.compound_id

            # Skip control compounds
            if parsed_result.is_control:
                self.log(f"  Skipping control compound: {compound_id}")
                continue

            # Try to match compound
            compound = self._match_compound(compound_id)

            if compound:
                matched_count += 1
                self.log(f"  Matched: {compound_id} -> {compound.reg_number}")
                results_to_import.append((parsed_result, compound))
            else:
                unmatched_count += 1
                if skip_unmatched:
                    self.log(f"  Skipping unmatched: {compound_id}")
                else:
                    self.stdout.write(self.style.WARNING(f"  Unmatched compound: {compound_id}"))
                    results_to_import.append((parsed_result, None))

        self.stdout.write(f"\nCompound matching:")
        self.stdout.write(f"  Matched: {matched_count}")
        self.stdout.write(f"  Unmatched: {unmatched_count}")

        # Import results
        if not preview and results_to_import:
            with transaction.atomic():
                self._import_results(
                    results_to_import,
                    parser=parser,
                    filepath=filepath,
                    metadata=parse_result.metadata,
                )
            self.stdout.write(self.style.SUCCESS(f"\nImported {len(results_to_import)} results"))
        elif preview:
            self._preview_results(results_to_import, parser)

        return {
            'results_count': len(parse_result.results),
            'errors_count': len([e for e in parse_result.errors if e.severity == 'error']),
            'matched_count': matched_count,
            'unmatched_count': unmatched_count,
        }

    def _match_compound(self, compound_id: str) -> Compound | None:
        """Try to match a compound ID to a registered compound.

        Supports:
        - NCL-XXXXXXXX format (direct reg_number match)
        - Plain numeric IDs (tries NCL-{padded})
        """
        if not compound_id:
            return None

        # Normalize compound ID
        compound_id = str(compound_id).strip()

        # Direct match on reg_number
        try:
            return Compound.objects.get(reg_number__iexact=compound_id)
        except Compound.DoesNotExist:
            pass

        # Try NCL format if numeric
        if compound_id.isdigit():
            ncl_id = f"NCL-{compound_id.zfill(8)}"
            try:
                return Compound.objects.get(reg_number__iexact=ncl_id)
            except Compound.DoesNotExist:
                pass

        # Try extracting number from NCL-XXXXX pattern
        match = re.match(r'NCL-?(\d+)', compound_id, re.IGNORECASE)
        if match:
            number = match.group(1)
            ncl_id = f"NCL-{number.zfill(8)}"
            try:
                return Compound.objects.get(reg_number__iexact=ncl_id)
            except Compound.DoesNotExist:
                pass

        return None

    def _get_or_create_protocol(self, parser) -> Protocol:
        """Get or create the protocol for this assay type."""
        protocol, created = Protocol.objects.get_or_create(
            name=f"{parser.vendor} {parser.assay_type}",
            defaults={
                'analysis_method': 'table_of_values',
                'comments': f"Auto-created for {parser.vendor} {parser.assay_type} imports",
            }
        )
        if created:
            self.stdout.write(f"Created protocol: {protocol.name}")
        return protocol

    def _import_results(
        self,
        results_to_import: list,
        parser,
        filepath: Path,
        metadata: dict,
    ) -> None:
        """Import parsed results into the database."""

        # Get or create protocol
        protocol = self._get_or_create_protocol(parser)

        # Create Assay for this import
        assay = Assay.objects.create(
            protocol=protocol,
            comments=f"Imported from {filepath.name}",
        )
        self.log(f"Created assay: {assay.id}")

        # Create DataSeries and AnalysisResult for each compound
        for idx, (parsed_result, compound) in enumerate(results_to_import):
            # Create AnalysisResult with parsed data
            analysis = AnalysisResult.objects.create(
                status='valid',
                results=parsed_result.results,
            )

            # Create DataSeries
            data_series = DataSeries.objects.create(
                assay=assay,
                compound=compound,
                compound_name=parsed_result.compound_id,
                row=idx,
                start_column=0,
                end_column=0,
                extracted_data=parsed_result.results.get('time_course', {}),
                analysis=analysis,
            )

            self.log(f"Created data series for {parsed_result.compound_id}")

    def _preview_results(self, results_to_import: list, parser) -> None:
        """Preview what would be imported."""
        self.stdout.write("\nPreview of results to import:")
        self.stdout.write("-" * 60)

        for parsed_result, compound in results_to_import[:10]:  # Show first 10
            compound_str = compound.reg_number if compound else "[UNMATCHED]"
            kpi_field = parser.kpi_field
            kpi_value = parsed_result.results.get(kpi_field, 'N/A')

            self.stdout.write(
                f"  {parsed_result.compound_id:<20} -> {compound_str:<15} "
                f"{kpi_field}={kpi_value}"
            )

        if len(results_to_import) > 10:
            self.stdout.write(f"  ... and {len(results_to_import) - 10} more")

    def log(self, message: str) -> None:
        """Log a message if verbose mode is enabled."""
        if self.verbose:
            self.stdout.write(f"  {message}")
