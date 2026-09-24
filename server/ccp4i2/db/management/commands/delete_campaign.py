"""
Delete a campaign by name: the group, its sites, its projects and their files.

Why this exists: a campaign is a group plus a parent project plus one project
per dataset, so undoing one by hand means deleting the group and then each
project in turn, remembering to ask for the files each time. That is tedious
for a real campaign and worse for the throwaway ones make_demo_campaign
builds, which are made to be rebuilt.

What goes:

* the ProjectGroup, and with it (by cascade) its memberships, its
  CampaignSites and their site evaluations;
* every project in the group -- parent and members -- with its jobs, files,
  file uses and imports, exactly as ``DELETE /projects/<id>/`` removes them;
* each project's directory on disk, unless --keep-files. Removal goes through
  remove_project_directory, which refuses anything that does not look like a
  CCP4i2 project directory, so a stale or mistaken path cannot take anything
  else with it.

What stays: a project that ALSO belongs to another group. Deleting it would
damage a campaign nobody asked to delete, so it only loses its membership of
this one, and the command says so.

Nothing is deleted without confirmation: the campaign name has to be typed
back, or --yes given. --dry-run reports the plan and touches nothing.

Usage:
    python manage.py delete_campaign BAZ2B_demo_campaign --dry-run
    python manage.py delete_campaign BAZ2B_demo_campaign
    python manage.py delete_campaign BAZ2B_demo_campaign --yes
    python manage.py delete_campaign BAZ2B_demo_campaign --keep-files
"""

import logging
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = (
        "Delete a campaign (project group) by name, with its sites, its "
        "parent and member projects, and their files on disk."
    )

    def add_arguments(self, parser):
        parser.add_argument("name", help="Campaign (project group) name.")
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Report what would be deleted; delete nothing.",
        )
        parser.add_argument(
            "--yes",
            action="store_true",
            help="Do not ask for confirmation.",
        )
        parser.add_argument(
            "--keep-files",
            action="store_true",
            help=(
                "Delete the database records only and leave the project "
                "directories on disk, from where they can be re-imported."
            ),
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="Delete even if some of the campaign's jobs are queued or running.",
        )

    def handle(self, *args, **options):
        from ccp4i2.db import models
        from ccp4i2.db.delete_project_directory import remove_project_directory

        name = options["name"]
        try:
            group = models.ProjectGroup.objects.get(name=name)
        except models.ProjectGroup.DoesNotExist:
            known = ", ".join(
                models.ProjectGroup.objects.order_by("name").values_list(
                    "name", flat=True
                )
            )
            raise CommandError(
                f"No campaign called '{name}'. "
                + (f"Campaigns: {known}" if known else "There are no campaigns.")
            )

        # Parent first, so the report reads in the order the campaign does.
        memberships = sorted(
            group.memberships.select_related("project"),
            key=lambda m: (
                m.type != models.ProjectGroupMembership.MembershipType.PARENT,
                m.project.name,
            ),
        )
        doomed, shared = [], []
        for membership in memberships:
            project = membership.project
            elsewhere = list(
                project.group_memberships.exclude(group=group).values_list(
                    "group__name", flat=True
                )
            )
            (shared if elsewhere else doomed).append((membership, elsewhere))

        keep_files = options["keep_files"]
        self.stdout.write(f"Campaign '{group.name}' ({group.get_type_display()})")
        self.stdout.write(
            f"  sites: {group.site_set.count()}, site evaluations: "
            f"{models.SiteEvaluation.objects.filter(site__group=group).count()}"
        )
        for membership, _ in doomed:
            project = membership.project
            self.stdout.write(
                f"  {membership.type:<6} {project.name}: "
                f"{project.jobs.count()} job(s), "
                f"{models.File.objects.filter(job__project=project).count()} "
                f"file(s), {project.directory}"
                + (" (directory kept)" if keep_files else "")
            )
        for membership, elsewhere in shared:
            self.stdout.write(
                self.style.WARNING(
                    f"  {membership.type:<6} {membership.project.name}: KEPT, "
                    f"it is also in {', '.join(elsewhere)}"
                )
            )

        active = models.Job.objects.filter(
            project__in=[m.project for m, _ in doomed],
            status__in=[
                models.Job.Status.QUEUED,
                models.Job.Status.RUNNING,
                models.Job.Status.RUNNING_REMOTELY,
            ],
        )
        if active.exists():
            listing = ", ".join(
                f"{job.project.name} job {job.number} ({job.get_status_display()})"
                for job in active.select_related("project")[:10]
            )
            if not options["force"]:
                raise CommandError(
                    f"Jobs are still active: {listing}. Deleting their "
                    "directories from under them leaves processes writing "
                    "into nothing. Wait, or use --force."
                )
            self.stdout.write(self.style.WARNING(f"  active jobs (--force): {listing}"))

        if options["dry_run"]:
            self.stdout.write(self.style.SUCCESS("\nDry run: nothing deleted."))
            return

        if not options["yes"]:
            typed = input(
                "\nThis cannot be undone. Type the campaign name to delete it: "
            )
            if typed.strip() != group.name:
                raise CommandError("Name did not match; nothing deleted.")

        directories = [Path(m.project.directory) for m, _ in doomed]
        with transaction.atomic():
            for membership, _ in doomed:
                project = membership.project
                # As ProjectViewSet.destroy: XData rows restrict the deletion
                # of their job; everything else hanging off a project cascades.
                models.XData.objects.filter(job__project=project).delete()
                project.delete()
            # Takes the memberships of any kept projects, the sites and their
            # evaluations with it.
            group.delete()
        logger.warning("Deleted campaign %s and %d project(s)", name, len(doomed))

        # After the records, and outside the transaction: a directory cannot be
        # rolled back, so none goes until the database side has committed.
        kept = []
        if not keep_files:
            for directory in directories:
                removed, reason = remove_project_directory(directory)
                if not removed:
                    kept.append((directory, reason))

        self.stdout.write(
            self.style.SUCCESS(
                f"\nDeleted campaign '{name}' and {len(doomed)} project(s)."
            )
        )
        for directory, reason in kept:
            self.stdout.write(
                self.style.WARNING(f"  left on disk: {directory} ({reason})")
            )
