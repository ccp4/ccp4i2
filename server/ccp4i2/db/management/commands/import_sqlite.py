"""
Django management command to import a legacy CCP4i2 SQLite database.

On --dry-run, also runs filesystem validation (project dirs, job dirs, files).

Usage:
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite --dry-run --verbose
    ccp4-python manage.py import_sqlite ~/.CCP4I2/db/database.sqlite --remap-dirs /old/path /new/path
"""

from django.core.management.base import BaseCommand, CommandError

from ccp4i2.db.import_sqlite import (
    SQLiteImporter,
    SQLiteValidator,
    StructuralIssuesError,
)


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
        parser.add_argument(
            "--copy-files",
            action="store_true",
            help="Copy each legacy project tree into the projects root and "
                 "repoint it (self-contained), instead of adopting it in place",
        )
        parser.add_argument(
            "--dest-root",
            default=None,
            help="Destination root for copied projects "
                 "(default: settings.CCP4I2_PROJECTS_DIR)",
        )
        parser.add_argument(
            "--allow-structural-issues",
            action="store_true",
            help="Proceed even if blocking structural issues (e.g. destination "
                 "collisions) are present",
        )

    def handle(self, *args, **options):
        remap_dirs = tuple(options["remap_dirs"]) if options["remap_dirs"] else None
        copy_files = options["copy_files"]
        dest_root = options["dest_root"]

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
                copy_files=copy_files,
                dest_root=dest_root,
            )
            report = validator.run()

            summary = report["summary"]
            self.stdout.write(f"\n  Projects on disk:       {summary['projects_on_disk']}")
            self.stdout.write(f"  Jobs on disk:           {summary['jobs_on_disk']}")
            self.stdout.write(f"  Files on disk:          {summary['files_on_disk']}")
            self.stdout.write(f"  Import sources on disk: {summary['import_sources_on_disk']}")
            self.stdout.write(f"  Integrity issues:       {summary['integrity_issues']}")
            self.stdout.write(f"  Data quality issues:    {summary['data_quality_issues']}")
            self.stdout.write(f"  Structure issues:       {summary['structure_issues']} "
                              f"({summary['blocking_issues']} blocking)")
            ps = summary["plan_summary"]
            self.stdout.write(f"  Plan: {ps['in_place']} in place, {ps['copy']} copied "
                              f"({ps['copied_due_to_nesting']} due to nesting)")

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
                copy_files=copy_files,
                dest_root=dest_root,
                allow_structural_issues=options["allow_structural_issues"],
            )

            try:
                result = importer.run()
            except StructuralIssuesError as e:
                self.stderr.write(self.style.ERROR(
                    "Refusing to import: blocking structural issues found."
                ))
                for issue in e.structure["issues"]:
                    if issue["severity"] == "blocking" and not issue.get("resolution"):
                        self.stderr.write(f"  - {issue['detail']}")
                raise CommandError(
                    "Re-run with --allow-structural-issues to proceed anyway, "
                    "or fix the issues above."
                )

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
