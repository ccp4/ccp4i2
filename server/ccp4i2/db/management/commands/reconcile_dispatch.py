"""Reconcile jobs whose program runs on a run target (RUNNING_REMOTELY).

    manage.py reconcile_dispatch --job <uuid>
    manage.py reconcile_dispatch --all [--json]

Idempotent: a job whose run is still going is left alone; one whose run has
ended is started again to complete from it. The same function the job page
calls (lib/utils/jobs/dispatch_record.reconcile).
"""
import json

from django.core.management.base import BaseCommand, CommandError

from ccp4i2.db import models
from ccp4i2.lib.utils.jobs.dispatch_record import reconcile, reconcile_all


class Command(BaseCommand):
    help = "Reconcile jobs running remotely on a program run target"

    def add_arguments(self, parser):
        parser.add_argument("--job", help="Job uuid")
        parser.add_argument("--all", action="store_true", help="Every job in RUNNING_REMOTELY")
        parser.add_argument("--json", action="store_true", help="One JSON object per line")

    def handle(self, *args, **options):
        if bool(options["job"]) == bool(options["all"]):
            raise CommandError("give --job <uuid> or --all")
        if options["job"]:
            try:
                job = models.Job.objects.get(uuid=options["job"])
            except (models.Job.DoesNotExist, ValueError):
                raise CommandError(f"no job with uuid {options['job']}")
            results = [(job, reconcile(job))]
        else:
            results = reconcile_all()
        for job, result in results:
            if options["json"]:
                self.stdout.write(json.dumps({"job": str(job.uuid), "number": job.number, **result}))
            else:
                self.stdout.write(f"job {job.number} ({job.uuid}): {result['action']} -- {result['reason']}")
