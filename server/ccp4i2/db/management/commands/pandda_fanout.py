"""Fan a PanDDA output tree out into per-dataset receipt jobs, from a shell.

The implementation is ``lib/utils/jobs/pandda_fanout.py``, shared with the
``pandda_fanout`` task; this command is for scripting and for a tree that
arrived by hand. Usage::

    manage.py pandda_fanout --tree <pandda2_out> --manifest <manifest.json> [--dry-run] [--no-run] [--json]
"""
import json
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError

from ccp4i2.lib.utils.jobs.pandda_fanout import (
    PROCESSED_DIR, Plan, execute_fanout, plan_fanout, read_manifest)


class Command(BaseCommand):
    help = "Fan a PanDDA output tree out into per-dataset pandda_events receipts"

    def add_arguments(self, parser):
        parser.add_argument("--tree", required=True, type=Path,
                            help="The pandda2_out directory (analyses/, processed_datasets/)")
        parser.add_argument("--manifest", required=True, type=Path,
                            help="manifest.json the orchestrator wrote when it staged the run")
        parser.add_argument("--run-job", default=None,
                            help="UUID of the pandda_campaign job; default: from the manifest; "
                                 "omit for a tree produced elsewhere")
        parser.add_argument("--dry-run", action="store_true", help="Report the plan; create nothing")
        parser.add_argument("--no-run", action="store_true",
                            help="Create the receipts but do not run them")
        parser.add_argument("--json", action="store_true", help="Report as JSON")

    def handle(self, **options):
        tree = options["tree"].expanduser().resolve()
        if not (tree / PROCESSED_DIR).is_dir():
            raise CommandError(f"{tree} has no {PROCESSED_DIR}/: not a PanDDA output tree")
        try:
            manifest = read_manifest(options["manifest"].expanduser())
        except (OSError, ValueError, json.JSONDecodeError) as err:
            raise CommandError(f"cannot read manifest: {err}")
        run_job_uuid = options["run_job"] or (manifest.get("provenance") or {}).get("run_job_uuid")

        plan = plan_fanout(tree, manifest, run_job_uuid)
        if not options["dry_run"]:
            execute_fanout(plan, run=not options["no_run"])
        self._report(plan, options)

    def _report(self, plan: Plan, options):
        if options["json"]:
            self.stdout.write(json.dumps({
                "tree": str(plan.tree), "run_job_uuid": plan.run_job_uuid,
                "incomplete": plan.incomplete, "dry_run": options["dry_run"],
                "outcomes": [o.as_dict() for o in plan.outcomes],
            }, indent=2))
            return
        mode = "DRY RUN: " if options["dry_run"] else ""
        self.stdout.write(f"{mode}{plan.tree}")
        self.stdout.write(f"run job: {plan.run_job_uuid or '(none: external tree)'}; "
                          f"events table: {'absent (partial run)' if plan.incomplete else 'present'}")
        for o in plan.outcomes:
            self.stdout.write(f"  {o.xtal:<12} {o.action:<10} {o.project or '-':<24} {o.label:<20} {o.reason}")
        summary = ", ".join(f"{plan.count(a)} {a}" for a in ("created", "skipped", "absent", "failed", "no_project")
                            if plan.count(a))
        self.stdout.write(summary or "nothing to do")
