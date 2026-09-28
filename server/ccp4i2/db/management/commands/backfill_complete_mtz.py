"""Register the complete reflection file of finished dimple jobs that ran
before the pipeline tracked it.

i2Dimple has registered dimple's unsplit ``final.mtz`` as its ``COMPLETE_MTZ``
output since 3.1.0a9 (PR #254, July 2026). Jobs run before then left the
file on disk in the i2Dimple job directory but registered only the model and
the two map-coefficient files, so anything that looks the reflections up by
their File row (the PanDDA campaign task's fill from campaign, the Export MTZ
button, purge protection) finds nothing and falls back to a bare path, which
the interface then shows as an empty field.

This command creates the missing row exactly as the gleaner would have:
against the i2Dimple job itself (never the SubstituteLigand pipeline that ran
it: the pipeline harvests its other outputs by copy but does not harvest this
one), in the job directory, so no file moves::

    job_param_name = COMPLETE_MTZ
    name           = final.mtz
    directory      = JOB_DIR (the i2Dimple job's own directory)
    type           = application/CCP4-mtz
    sub_type       = 0, content = 0   (what CMtzDataFile registers unset)
    annotation     = "Complete unsplit reflection file from dimple"

Dry run by default; nothing is written without ``--commit``.

Usage:
    python manage.py backfill_complete_mtz                     # dry run, every project
    python manage.py backfill_complete_mtz --commit
    python manage.py backfill_complete_mtz --campaign <campaign-name> --commit
    python manage.py backfill_complete_mtz --projectname <project-name>
    python manage.py backfill_complete_mtz --commit --limit 10  # a first, small batch
"""
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from ccp4i2.db.models import File, FileType, Job, Project, ProjectGroup, ProjectGroupMembership
from ccp4i2.lib.pandda_export import DIMPLE_TASK_NAMES

PARAM_NAME = "COMPLETE_MTZ"
FILE_NAME = "final.mtz"
FILE_TYPE = "application/CCP4-mtz"
ANNOTATION = "Complete unsplit reflection file from dimple"


class Command(BaseCommand):
    help = "Register final.mtz as COMPLETE_MTZ on finished dimple jobs that predate its tracking"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("--commit", action="store_true",
                            help="Write the rows (default is a dry run that only reports)")
        parser.add_argument("-pn", "--projectname", type=str,
                            help="Limit to one project")
        parser.add_argument("--campaign", type=str,
                            help="Limit to the member projects of one campaign (project group name)")
        parser.add_argument("--limit", type=int,
                            help="Stop after registering this many rows (a first batch)")
        parser.add_argument("--quiet", action="store_true",
                            help="Print the summary only, not one line per job")

    def handle(self, *args, **options):
        commit = options["commit"]
        if not commit:
            self.stdout.write(self.style.WARNING("DRY RUN: pass --commit to write to the database\n"))

        jobs = (Job.objects.filter(status=Job.Status.FINISHED, task_name__in=list(DIMPLE_TASK_NAMES))
                .select_related("project").order_by("project__name", "number"))
        if options["projectname"]:
            if not Project.objects.filter(name=options["projectname"]).exists():
                raise CommandError(f"Project '{options['projectname']}' not found")
            jobs = jobs.filter(project__name=options["projectname"])
        if options["campaign"]:
            try:
                group = ProjectGroup.objects.get(name=options["campaign"])
            except ProjectGroup.DoesNotExist:
                raise CommandError(f"Campaign '{options['campaign']}' not found")
            members = ProjectGroupMembership.objects.filter(
                group=group, type=ProjectGroupMembership.MembershipType.MEMBER
            ).values_list("project_id", flat=True)
            jobs = jobs.filter(project_id__in=list(members))

        already = set(File.objects.filter(job__in=jobs, job_param_name=PARAM_NAME)
                      .values_list("job_id", flat=True))
        to_register = []
        skipped_registered = 0
        skipped_no_file = 0
        for job in jobs.iterator():
            if job.id in already:
                skipped_registered += 1
                continue
            if not (job.directory / FILE_NAME).is_file():
                skipped_no_file += 1
                if not options["quiet"]:
                    self.stdout.write(f"  {job.project.name}/{job.number}: no {FILE_NAME} on disk, skipped")
                continue
            to_register.append(job)
            if not options["quiet"]:
                self.stdout.write(f"  {job.project.name}/{job.number}: {job.directory / FILE_NAME}")
            if options["limit"] and len(to_register) >= options["limit"]:
                break

        if commit and to_register:
            with transaction.atomic():
                file_type, _ = FileType.objects.get_or_create(
                    name=FILE_TYPE, defaults={"description": "MTZ experimental data"})
                File.objects.bulk_create([
                    File(name=FILE_NAME, directory=File.Directory.JOB_DIR, type=file_type,
                         sub_type=0, content=0, annotation=ANNOTATION,
                         job=job, job_param_name=PARAM_NAME)
                    for job in to_register
                ])

        verb = "Registered" if commit else "Would register"
        self.stdout.write("")
        self.stdout.write(self.style.SUCCESS(f"{verb} {len(to_register)} {PARAM_NAME} row(s)"))
        self.stdout.write(f"  Skipped (already registered): {skipped_registered}")
        self.stdout.write(f"  Skipped (no {FILE_NAME} on disk): {skipped_no_file}")
