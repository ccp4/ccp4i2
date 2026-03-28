# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Django management command to import a legacy CCP4i2 SQLite database.

On --dry-run, also runs filesystem validation (project dirs, job dirs, files).

Usage:
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite --dry-run --verbose
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite --remap-dirs /old/path /new/path
"""

from django.core.management.base import BaseCommand

from ccp4i2.db.import_sqlite import SQLiteImporter, SQLiteValidator


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
            help="Validate and simulate import without committing changes",
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
        remap_dirs = tuple(options["remap_dirs"]) if options["remap_dirs"] else None

        def log_fn(msg):
            self.stdout.write(msg)

        if options["dry_run"]:
            # Dry-run: validate only, no Django DB needed
            self.stdout.write(f"\nValidating: {options['db_path']}")
            self.stdout.write("=" * 60)

            validator = SQLiteValidator(
                db_path=options["db_path"],
                remap_dirs=remap_dirs,
                verbose=options["verbose"],
                log_fn=log_fn,
            )
            report = validator.run()

            summary = report["summary"]
            self.stdout.write(f"\n  Projects on disk:       {summary['projects_on_disk']}")
            self.stdout.write(f"  Jobs on disk:           {summary['jobs_on_disk']}")
            self.stdout.write(f"  Files on disk:          {summary['files_on_disk']}")
            self.stdout.write(f"  Import sources on disk: {summary['import_sources_on_disk']}")
            self.stdout.write(f"  Integrity issues:       {summary['integrity_issues']}")
            self.stdout.write(f"  Data quality issues:    {summary['data_quality_issues']}")

            self.stdout.write("\n" + "=" * 60)
            if summary["ok"]:
                self.stdout.write(self.style.SUCCESS("All validation checks passed"))
            else:
                self.stdout.write(self.style.WARNING("Some validation checks failed — review details above"))
            self.stdout.write("=" * 60)
        else:
            # Real import
            importer = SQLiteImporter(
                db_path=options["db_path"],
                remap_dirs=remap_dirs,
                dry_run=False,
                continue_on_error=options["continue_on_error"],
                verbose=options["verbose"],
                log_fn=log_fn,
            )

            result = importer.run()

            self.stdout.write("\n" + "-" * 60)
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
