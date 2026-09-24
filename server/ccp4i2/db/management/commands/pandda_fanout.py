"""Fan a PanDDA output tree out into per-dataset receipt jobs.

Design note: docs/pandda-campaign-design.md, section 8.

Fan-out takes a conformant ``pandda2_out/`` tree and the manifest the
orchestrator wrote when it staged the run, and creates one ``pandda_events``
receipt job in each dataset's own project. It does not care who produced
the tree: this machine, a cluster, or a colleague's disk (8.4). It calls the
same entry points the interface calls (create the job, set its parameters,
run it), never a parallel implementation of job creation (8.1).

Four properties, each the answer to a way the previous attempt failed (8.2):

* idempotent: keyed on (run job uuid, dtag), read from each existing receipt's
  own parameters, so running it again creates nothing that already landed;
* previewable: ``--dry-run`` reports the plan and touches nothing;
* reported: one row per dataset, created / skipped / absent / failed, with
  the reason;
* retryable: a failure is recorded against its dataset and the rest go on.

A partial tree is fine: datasets the run got to become receipts, the rest
are reported absent, and every receipt from a run without an events table
is told so (``RUN_INCOMPLETE``).

Usage::

    manage.py pandda_fanout --tree <pandda2_out> --manifest <manifest.json> [--dry-run] [--no-run]
"""
import json
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from django.core.management.base import BaseCommand, CommandError

from ccp4i2.db import models
from ccp4i2.wrappers.pandda_campaign.script.pandda_staging import read_manifest
from ccp4i2.wrappers.pandda_events.script.pandda_tree import (
    ANALYSES_DIR, EVENTS_TABLE, PROCESSED_DIR)

RECEIPT_TASK = "pandda_events"


@dataclass
class Outcome:
    xtal: str
    label: str
    action: str          # created | skipped | absent | failed | no_project
    reason: str = ""
    job_uuid: Optional[str] = None
    project: Optional[str] = None

    def as_dict(self):
        return {"xtal": self.xtal, "label": self.label, "action": self.action,
                "reason": self.reason, "job_uuid": self.job_uuid, "project": self.project}


@dataclass
class Plan:
    tree: Path
    run_job_uuid: Optional[str]
    incomplete: bool
    outcomes: List[Outcome] = field(default_factory=list)

    def count(self, action):
        return sum(1 for o in self.outcomes if o.action == action)


# ---------------------------------------------------------------------------
# Reading what exists
# ---------------------------------------------------------------------------

def tree_is_complete(tree: Path) -> bool:
    return (tree / ANALYSES_DIR / EVENTS_TABLE).is_file()


def _receipt_key(job: models.Job):
    """``(run_job_uuid, dtag, tree)`` as a receipt's own parameters record
    them, or None if they cannot be read."""
    for name in ("params.xml", "input_params.xml"):
        path = Path(job.directory) / name
        if not path.is_file():
            continue
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            continue
        inputs = root.find("ccp4i2_body/inputData")
        if inputs is None:
            continue
        get = lambda tag: (inputs.findtext(tag) or "").strip()
        return get("RUN_JOB_UUID"), get("DTAG"), get("PANDDA_OUT_DIR")
    return None


def existing_receipt(project: models.Project, run_job_uuid: Optional[str],
                     dtag: str, tree: Path) -> Optional[models.Job]:
    """The receipt already in ``project`` for this (run, dtag), if any.

    With a run job the key is (run job uuid, dtag). A tree produced elsewhere
    has no run job, so the key falls back to (tree path, dtag).
    """
    for job in models.Job.objects.filter(project=project, task_name=RECEIPT_TASK).exclude(
            status=models.Job.Status.TO_DELETE).order_by("id"):
        key = _receipt_key(job)
        if key is None:
            continue
        run, job_dtag, job_tree = key
        if job_dtag != dtag:
            continue
        if run_job_uuid and run.replace("-", "") == run_job_uuid.replace("-", ""):
            return job
        if not run_job_uuid and not run and Path(job_tree or "") == tree:
            return job
    return None


# ---------------------------------------------------------------------------
# Planning and doing
# ---------------------------------------------------------------------------

def plan_fanout(tree: Path, manifest: dict, run_job_uuid: Optional[str]) -> Plan:
    """What fan-out would do, dataset by dataset, without doing it."""
    tree = Path(tree)
    plan = Plan(tree=tree, run_job_uuid=run_job_uuid, incomplete=not tree_is_complete(tree))
    for entry in manifest.get("datasets", []):
        xtal, label = entry["xtal"], entry.get("label", "")
        project_uuid = entry.get("project_uuid")
        project = models.Project.objects.filter(uuid=project_uuid).first() if project_uuid else None
        if project is None:
            plan.outcomes.append(Outcome(xtal, label, "no_project",
                                         f"project {project_uuid or '(none recorded)'} is not in this database"))
            continue
        if not (tree / PROCESSED_DIR / xtal).is_dir():
            plan.outcomes.append(Outcome(xtal, label, "absent",
                                         "the run wrote nothing for this dataset", project=project.name))
            continue
        existing = existing_receipt(project, run_job_uuid, xtal, tree)
        if existing is not None:
            plan.outcomes.append(Outcome(xtal, label, "skipped",
                                         f"receipt already exists: job {existing.number}",
                                         job_uuid=str(existing.uuid), project=project.name))
            continue
        plan.outcomes.append(Outcome(xtal, label, "created", "", project=project.name))
    return plan


def _run_label(run_job_uuid: Optional[str]) -> str:
    if not run_job_uuid:
        return "an external PanDDA run"
    job = models.Job.objects.filter(uuid=run_job_uuid).select_related("project").first()
    if job is None:
        return f"PanDDA run {run_job_uuid[:8]}"
    return f"PanDDA run {job.project.name}/{job.number}"


def create_receipt(project: models.Project, xtal: str, label: str, tree: Path,
                   run_job_uuid: Optional[str], incomplete: bool) -> models.Job:
    """One receipt job, parameterised, through the same helpers the interface
    uses. Raises on any parameter that could not be set."""
    from ccp4i2.lib.utils.jobs.create import create_job
    from ccp4i2.lib.utils.parameters.set_param import set_parameter

    title = f"PanDDA events for {label or xtal} from {_run_label(run_job_uuid)}"
    job_uuid = create_job(projectId=str(project.uuid), taskName=RECEIPT_TASK, title=title[:255])
    job = models.Job.objects.get(uuid=job_uuid)
    values = [
        ("inputData.PANDDA_OUT_DIR", str(tree)),
        ("inputData.DTAG", xtal),
        ("inputData.RUN_INCOMPLETE", bool(incomplete)),
    ]
    if run_job_uuid:
        values.append(("inputData.RUN_JOB_UUID", str(run_job_uuid)))
    for path, value in values:
        result = set_parameter(job, path, value)
        if not getattr(result, "success", False):
            raise CommandError(f"could not set {path}={value!r} on job {job.number}: "
                               f"{getattr(result, 'error', 'unknown error')}")
    return job


def execute_fanout(plan: Plan, run: bool = True, stdout=None) -> Plan:
    """Do what the plan says. A failure on one dataset is recorded and the
    rest continue; nothing is rolled back, so a rerun picks up the remainder."""
    for outcome in plan.outcomes:
        if outcome.action != "created":
            continue
        project = models.Project.objects.get(name=outcome.project)
        try:
            job = create_receipt(project, outcome.xtal, outcome.label, plan.tree,
                                 plan.run_job_uuid, plan.incomplete)
        except Exception as err:      # noqa: BLE001 - recorded, not raised
            outcome.action = "failed"
            outcome.reason = f"could not create the receipt: {err}"
            continue
        outcome.job_uuid = str(job.uuid)
        outcome.reason = f"job {job.number}"
        if not run:
            continue
        from ccp4i2.lib.utils.jobs.context_run import run_job_context_aware
        result = run_job_context_aware(job, synchronous=True)
        if not result.get("success"):
            outcome.action = "failed"
            outcome.reason = f"job {job.number} created but did not start: {result.get('error')}"
            continue
        job.refresh_from_db()
        outcome.reason = f"job {job.number}: {job.get_status_display()}"
    return plan


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
            execute_fanout(plan, run=not options["no_run"], stdout=self.stdout)
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
