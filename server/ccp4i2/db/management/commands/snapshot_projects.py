"""Write the on-disk recovery snapshot for projects that have none.

Projects created before snapshots existed never trigger one: a finished,
dormant project is both the kind you most want protected and the kind nothing
will touch again. This retrofits them.

Cheap enough to run on every start-up in its default ``--missing-only`` mode,
which skips projects already covered.
"""

from django.core.management.base import BaseCommand

from ...project_snapshot import snapshot_all_projects, snapshot_status


class Command(BaseCommand):
    help = "Write DATABASE.db.xml recovery snapshots for projects on disk."

    def add_arguments(self, parser):
        parser.add_argument(
            "--all",
            action="store_true",
            help="Rewrite every project's snapshot, not just the missing ones.",
        )
        parser.add_argument(
            "--status",
            action="store_true",
            help="Report which projects are protected, and write nothing.",
        )

    def handle(self, *args, **options):
        if options["status"]:
            self._report_status()
            return

        result = snapshot_all_projects(missing_only=not options["all"])

        self.stdout.write(
            f"Wrote {result['written']} snapshot(s); "
            f"{result['already_current']} already current, "
            f"{result['skipped']} could not be written."
        )

    def _report_status(self):
        rows = snapshot_status()
        if not rows:
            self.stdout.write("No projects.")
            return

        width = max(len(row["name"]) for row in rows)
        for row in rows:
            if not row["directory_exists"]:
                state = "directory missing"
            elif row["has_snapshot"]:
                state = "protected"
            elif row["foreign_snapshot"]:
                state = "has a snapshot from another installation"
            else:
                state = "NOT PROTECTED"
            self.stdout.write(f"  {row['name']:<{width}}  {state}")

        unprotected = sum(
            1 for row in rows if row["directory_exists"] and not row["has_snapshot"]
        )
        self.stdout.write("")
        self.stdout.write(
            f"{len(rows) - unprotected} of {len(rows)} project(s) protected."
        )
        if unprotected:
            self.stdout.write("Run without --status to write the missing snapshots.")
