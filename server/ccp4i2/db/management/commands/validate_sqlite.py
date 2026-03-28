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
Django management command to validate a legacy CCP4i2 SQLite database.

Checks that project directories, job directories, files, and imported file
sources all exist on disk. Also checks referential integrity and data quality.

Usage:
    ccp4-python manage.py validate_sqlite ~/.CCP4I2/db/database.sqlite
    ccp4-python manage.py validate_sqlite ~/.CCP4I2/db/database.sqlite --verbose
    ccp4-python manage.py validate_sqlite ~/.CCP4I2/db/database.sqlite --remap-dirs /old /new
    ccp4-python manage.py validate_sqlite ~/.CCP4I2/db/database.sqlite --json
"""

import json

from django.core.management.base import BaseCommand

from ccp4i2.db.import_sqlite import SQLiteValidator


class Command(BaseCommand):
    help = "Validate a legacy CCP4i2 SQLite database against the filesystem"

    def add_arguments(self, parser):
        parser.add_argument(
            "db_path",
            type=str,
            help="Path to legacy CCP4i2 SQLite database (e.g. ~/.CCP4I2/db/database.sqlite)",
        )
        parser.add_argument(
            "--verbose",
            action="store_true",
            help="List individual missing paths",
        )
        parser.add_argument(
            "--remap-dirs",
            nargs=2,
            metavar=("FROM", "TO"),
            help="Remap project directories before checking (e.g., --remap-dirs /old/path /new/path)",
        )
        parser.add_argument(
            "--json",
            action="store_true",
            dest="output_json",
            help="Output full report as JSON",
        )

    def handle(self, *args, **options):
        validator = SQLiteValidator(
            db_path=options["db_path"],
            remap_dirs=tuple(options["remap_dirs"]) if options["remap_dirs"] else None,
            verbose=options["verbose"],
            log_fn=lambda msg: self.stdout.write(msg),
        )

        self.stdout.write(f"\nValidating: {options['db_path']}")
        self.stdout.write("-" * 60)

        report = validator.run()

        if options["output_json"]:
            self.stdout.write(json.dumps(report, indent=2))
            return

        self.stdout.write("\n" + "-" * 60)
        self.stdout.write("Summary")
        self.stdout.write("-" * 60)

        summary = report["summary"]
        self.stdout.write(f"  Projects on disk:       {summary['projects_on_disk']}")
        self.stdout.write(f"  Jobs on disk:           {summary['jobs_on_disk']}")
        self.stdout.write(f"  Files on disk:          {summary['files_on_disk']}")
        self.stdout.write(f"  Import sources on disk: {summary['import_sources_on_disk']}")
        self.stdout.write(f"  Integrity issues:       {summary['integrity_issues']}")
        self.stdout.write(f"  Data quality issues:    {summary['data_quality_issues']}")

        if summary["ok"]:
            self.stdout.write(self.style.SUCCESS("\n  All checks passed"))
        else:
            self.stdout.write(self.style.WARNING("\n  Some checks failed — review details above"))

        self.stdout.write("-" * 60)
