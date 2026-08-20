"""Rebuild the database from the snapshots left in project directories.

Restore, not import: the directories are the truth and the database is being
repaired, so nothing is renumbered and every job must come back under the
number its directory already has on disk.

    manage.py restore_projects --registry            # everything the registry knows
    manage.py restore_projects --scan ~/CCP4I2_PROJECTS
    manage.py restore_projects --directory ~/CCP4I2_PROJECTS/BetaBlip
    manage.py restore_projects --registry --dry-run  # look first

An existing project is left alone unless --replace is given: overwriting a live
project with an older snapshot would lose work rather than recover it.
"""

from pathlib import Path

from django.core.management.base import BaseCommand, CommandError

from ...project_snapshot import registry_path
from ...restore_project import (
    discover_restorable,
    inspect,
    restore_all,
    restore_from_registry,
)


class Command(BaseCommand):
    help = "Rebuild database rows from the DATABASE.db.xml in project directories."

    def add_arguments(self, parser):
        source = parser.add_mutually_exclusive_group(required=True)
        source.add_argument(
            "--registry",
            action="store_true",
            help="Restore every project listed in the project directory registry.",
        )
        source.add_argument(
            "--scan",
            metavar="DIRECTORY",
            help="Restore every project directory found under DIRECTORY. "
            "For when the registry has been lost as well.",
        )
        source.add_argument(
            "--directory",
            metavar="DIRECTORY",
            help="Restore this one project directory.",
        )
        parser.add_argument(
            "--replace",
            action="store_true",
            help="Rebuild projects already in the database, discarding the rows "
            "currently held for them.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Report what would be restored and change nothing.",
        )

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        replace = options["replace"]

        if options["registry"]:
            self.stdout.write(f"Reading {registry_path()}")
            result = restore_from_registry(replace=replace, dry_run=dry_run)
            if result["registry_entries"] == 0:
                raise CommandError(
                    "The registry is empty or missing. Use --scan to look for "
                    "project directories directly."
                )
        elif options["scan"]:
            root = Path(options["scan"]).expanduser()
            directories = discover_restorable(root)
            if not directories:
                raise CommandError(f"No project directories with a snapshot under {root}")
            self.stdout.write(f"Found {len(directories)} project director(y/ies) under {root}")
            result = restore_all(directories, replace=replace, dry_run=dry_run)
        else:
            directory = Path(options["directory"]).expanduser()
            report = inspect(directory)
            if report.skipped_reason and report.project_uuid is None:
                raise CommandError(f"{directory}: {report.skipped_reason}")
            result = restore_all([directory], replace=replace, dry_run=dry_run)

        self._report(result, dry_run)

    def _report(self, result, dry_run):
        restored = result["restored"]
        skipped = result["skipped"]

        if dry_run:
            self.stdout.write("")
            self.stdout.write("Would restore:")
            for entry in skipped + restored:
                if entry["skipped_reason"] and entry["project_uuid"] is None:
                    continue
                if entry["skipped_reason"]:
                    continue
                self._describe(entry, prefix="  ")
            self.stdout.write("")

        for entry in restored:
            self._describe(entry, prefix="  restored ")

        for entry in skipped:
            if dry_run and not entry["skipped_reason"]:
                continue
            name = entry["project_name"] or entry["directory"]
            self.stdout.write(f"  skipped   {name}: {entry['skipped_reason']}")

        self.stdout.write("")
        verb = "would be restored" if dry_run else "restored"
        count = len([e for e in skipped + restored if not e["skipped_reason"]])
        self.stdout.write(
            f"{count if dry_run else len(restored)} project(s) {verb}, "
            f"{len([e for e in skipped if e['skipped_reason']])} skipped."
        )

    def _describe(self, entry, prefix=""):
        line = (
            f"{prefix}{entry['project_name']}: "
            f"{entry['jobs']} job(s), {entry['files']} file(s)"
        )
        if entry["relocated"]:
            line += f" [moved from {entry['recorded_directory']}]"
        self.stdout.write(line)
        for warning in entry["warnings"]:
            self.stdout.write(f"{' ' * len(prefix)}  warning: {warning}")
        if entry["missing_job_directories"]:
            numbers = ", ".join(entry["missing_job_directories"][:8])
            self.stdout.write(
                f"{' ' * len(prefix)}  job directories not on disk: {numbers}"
            )
