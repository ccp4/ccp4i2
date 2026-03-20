"""
Django management command to import a legacy CCP4i2 SQLite database.

Usage:
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite --dry-run --verbose
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite --remap-dirs /old/path /new/path
"""

from django.core.management.base import BaseCommand

from ccp4i2.db.import_sqlite import SQLiteImporter


class Command(BaseCommand):
    help = "Import legacy CCP4i2 data directly from a SQLite database file"

    def add_arguments(self, parser):
        parser.add_argument(
            "db_path",
            type=str,
            help="Path to legacy CCP4i2 SQLite database (e.g. ~/.CCP4I2/db/database.sqlite)",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Validate without committing changes",
        )
        parser.add_argument(
            "--verbose",
            action="store_true",
            help="Show detailed progress",
        )
        parser.add_argument(
            "--remap-dirs",
            nargs=2,
            metavar=("FROM", "TO"),
            help="Remap project directories (e.g., --remap-dirs /old/path /new/path)",
        )
        parser.add_argument(
            "--continue-on-error",
            action="store_true",
            help="Continue importing remaining records if one fails",
        )

    def handle(self, *args, **options):
        importer = SQLiteImporter(
            db_path=options["db_path"],
            remap_dirs=tuple(options["remap_dirs"]) if options["remap_dirs"] else None,
            dry_run=options["dry_run"],
            continue_on_error=options["continue_on_error"],
            verbose=options["verbose"],
            log_fn=lambda msg: self.stdout.write(msg),
        )

        result = importer.run()

        self.stdout.write("\n" + "-" * 60)
        if result["dry_run"]:
            self.stdout.write(self.style.WARNING("DRY RUN completed (no changes saved)"))
        else:
            self.stdout.write(self.style.SUCCESS("Import completed successfully!"))

        stats = result["stats"]
        self.stdout.write(f"  Total records processed: {sum(stats.values())}")
        for key, val in sorted(stats.items()):
            if val:
                self.stdout.write(f"    {key}: {val}")

        if result["errors"]:
            self.stdout.write(self.style.ERROR(f"\n  Errors: {len(result['errors'])}"))
            for err in result["errors"][:20]:
                self.stderr.write(f"    {err}")

        self.stdout.write("-" * 60)
