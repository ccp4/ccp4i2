"""Give a FileUse row that records a bare list index the list's name back.

A ``CList`` names its elements ``"[i]"``, and until PR #668 the import path
recorded that bare in ``FileUse.job_param_name``. So which list a file went into
was lost, and worse, one ``(job, job_param_name)`` pair came to name several
different files: on one real job, ``'[0]'`` named a PDB, a map and a difference
map at once. Measured on a live database: 162 such rows across 30 jobs, 44
colliding groups.

#668 stops new ones being written. This repairs the ones already there.

Deliberately a management command and NOT a data migration. The desktop runs
``migrate`` on every launch (client/main/ccp4i2-django-server.ts), so a repair
that trips over old data would fail the launch path, on exactly the installs most
likely to have odd data. This also needs to import a task's plugin to read its
parameter declarations, which for some tasks needs CCP4 or libtbx, and a
migration must never depend on that.

**Nothing is changed unless the answer is certain.** Ambiguity, missing
provenance, or a params file that disagrees all mean the row is left exactly as
it is and counted. A wrong provenance label is worse than a bare one.

Scoped by the broken SHAPE -- ``job_param_name`` matching ``^\\[\\d+\\]$`` -- not
by task, age or project. Qt-era ccp4i2 had no separate name for a CList member at
all, so rows imported from those databases carry the plain list name repeated;
they do not match this filter and are never touched. That is deliberate: they
resolve correctly already (the reference syntax falls back from
``XYZIN_LIST[1]`` to the name plus an index over the ordered files), and
"repairing" them would mean inventing indices.

A repaired name is still POSITIONAL, and therefore still not an identity key.
``upload_param.py`` matches a slot's current file by file identity rather than by
this name for exactly that reason -- a list insert or removal moves every index
after it. Do not start keying on ``(job, job_param_name)`` because it now looks
trustworthy.

Usage::

    manage.py repair_file_use_names                # dry run, the default
    manage.py repair_file_use_names --apply
    manage.py repair_file_use_names --project NAME --verbose
"""

import re
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from django.core.management.base import BaseCommand

BARE_INDEX = re.compile(r"^\[(\d+)\]$")


class Command(BaseCommand):
    help = (
        "Repair FileUse rows whose job_param_name is a bare list index, "
        "giving them the list's name. Dry run unless --apply."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--apply",
            action="store_true",
            help="Write the repairs. Without this nothing is changed.",
        )
        parser.add_argument(
            "--project",
            default=None,
            help="Restrict to one project, by name.",
        )
        parser.add_argument(
            "--verbose",
            action="store_true",
            help="Print a line per row, repaired or skipped.",
        )

    def handle(self, *args, **options):
        from ccp4i2.db import models

        apply_changes = options["apply"]
        verbose = options["verbose"]

        rows = models.FileUse.objects.filter(
            job_param_name__regex=r"^\[[0-9]+\]$"
        ).select_related("job", "job__project", "file", "file__type")
        if options["project"]:
            rows = rows.filter(job__project__name=options["project"])

        counts = Counter()
        reasons = Counter()
        self._lists_by_task = {}
        self._params_by_job = {}

        for row in rows.order_by("job_id", "id"):
            try:
                outcome, detail = self._consider(row)
            except Exception as err:  # noqa: BLE001 - never abort the run
                outcome, detail = "error", f"{type(err).__name__}: {err}"

            counts[outcome] += 1
            if outcome != "repair":
                reasons[detail] += 1
                if verbose:
                    self.stdout.write(
                        f"  skip  job {row.job.number} {row.job_param_name} "
                        f"{row.file.name}: {detail}"
                    )
                continue

            if verbose:
                self.stdout.write(
                    f"  {'repair' if apply_changes else 'would'}  job "
                    f"{row.job.number} {row.job_param_name} -> {detail} "
                    f"({row.file.name})"
                )
            if apply_changes:
                row.job_param_name = detail
                row.save(update_fields=["job_param_name"])

        self._report(counts, reasons, apply_changes)

    # ------------------------------------------------------------------
    def _consider(self, row):
        """Return ``("repair", new_name)`` or ``("skip", reason)``.

        Every path that cannot be certain returns a skip. The row is only
        repaired when the task declares exactly one list the file could belong
        to AND, where the job's params file can be read, that file really is in
        that list.
        """
        from ccp4i2.db import models

        index = BARE_INDEX.match(row.job_param_name).group(1)

        # The row's ROLE decides which half of the container can possibly hold
        # it: an IN use is an input list, an OUT use an output one. Without this
        # a pdb file was ambiguous between coot1's XYZIN_LIST and its XYZOUT
        # list, which left 66 of 162 rows unrepairable for no good reason.
        section = (
            "inputData" if row.role == models.FileUse.Role.IN else "outputData"
        )
        lists = self._lists_for(row.job.task_name, section)
        if lists is None:
            return "skip", f"cannot read the declarations of {row.job.task_name}"
        if not lists:
            return "skip", (
                f"{row.job.task_name} declares no file lists in {section}"
            )

        file_type = row.file.type.name  # type is non-nullable, unlike sub_type
        sub_type = str(row.file.sub_type) if row.file.sub_type is not None else None
        candidates = [
            name
            for name, mime, required_sub in lists
            if mime == file_type
            and (required_sub is None or required_sub == sub_type)
        ]
        if len(candidates) != 1:
            return "skip", (
                f"{len(candidates)} candidate lists in {section} for "
                f"{file_type} subType={sub_type} on {row.job.task_name}"
            )
        new_name = f"{candidates[0]}[{index}]"

        # Provenance: a file that is not on disk any more cannot have its
        # provenance confirmed, so its label is left alone.
        try:
            path = row.file.path
        except Exception:
            return "skip", "file path cannot be resolved"
        if not path or not Path(path).exists():
            return "skip", "file is missing from disk"

        # Where the job's params can be read they are the authority: the file
        # must actually appear in the list we are about to name. A params file
        # that cannot be read is not a reason to refuse, but a params file that
        # DISAGREES is.
        members = self._params_lists_for(row.job)
        if members is not None:
            db_file_id = str(row.file.uuid).replace("-", "")
            holding = [name for name, ids in members.items() if db_file_id in ids]
            if holding and candidates[0] not in holding:
                return "skip", (
                    f"params.xml puts this file in {holding}, not {candidates[0]}"
                )
            if not holding:
                return "skip", "params.xml does not list this file at all"

        if (
            models.FileUse.objects.filter(
                file=row.file,
                job=row.job,
                role=row.role,
                job_param_name=new_name,
            )
            .exclude(pk=row.pk)
            .exists()
        ):
            return "skip", f"a row for {new_name} already exists"

        return "repair", new_name

    # ------------------------------------------------------------------
    def _lists_for(self, task_name, section):
        """``[(list name, mimeType, requiredSubType)]`` for one container section.

        None means the task's declarations could not be read -- an unregistered
        task, or a plugin needing CCP4/libtbx that is not present. Cached,
        because importing a plugin is not cheap and jobs repeat.
        """
        key = (task_name, section)
        if key in self._lists_by_task:
            return self._lists_by_task[key]

        result = None
        try:
            import tempfile

            from ccp4i2.core.base_object.cdata_file import CDataFile
            from ccp4i2.core.base_object.fundamental_types import CList
            from ccp4i2.core.CCP4Container import CContainer
            from ccp4i2.core.tasks import get_plugin_class

            plugin = get_plugin_class(task_name)(
                workDirectory=tempfile.mkdtemp(), parent=None
            )
            found = []

            def walk(node):
                for child in node.children():
                    if isinstance(child, CList):
                        item = child.makeItem()
                        if isinstance(item, CDataFile):
                            found.append(
                                (
                                    child.objectName(),
                                    _qualifier(item, "mimeTypeName"),
                                    _qualifier(item, "requiredSubType"),
                                )
                            )
                    elif isinstance(child, CContainer):
                        walk(child)

            node = getattr(plugin.container, section, None)
            if node is not None:
                walk(node)
            result = found
        except Exception:
            result = None

        self._lists_by_task[key] = result
        return result

    def _params_lists_for(self, job):
        """``{list name: [dbFileId, ...]}`` from the job's params, or None.

        None means no params file could be read or parsed -- a relocated
        project, a legacy import, a job whose directory has gone. Not a reason
        to refuse a repair, only to do it unverified.
        """
        if job.id in self._params_by_job:
            return self._params_by_job[job.id]

        result = None
        try:
            directory = Path(job.directory)
            for name in ("params.xml", "input_params.xml"):
                candidate = directory / name
                if not candidate.exists():
                    continue
                root = ET.parse(candidate).getroot()
                members = {}
                for element in root.iter():
                    ids = [
                        node.text.strip().replace("-", "")
                        for node in element.iter("dbFileId")
                        if node.text and node.text.strip()
                    ]
                    if ids:
                        members.setdefault(element.tag, []).extend(ids)
                if members:
                    result = members
                    break
        except Exception:
            result = None

        self._params_by_job[job.id] = result
        return result

    # ------------------------------------------------------------------
    def _report(self, counts, reasons, apply_changes):
        total = sum(counts.values())
        self.stdout.write("")
        if not total:
            self.stdout.write("No FileUse row records a bare list index. Nothing to do.")
            return

        verb = "Repaired" if apply_changes else "Would repair"
        self.stdout.write(f"{total} row(s) with a bare list index:")
        self.stdout.write(f"  {verb}: {counts['repair']}")
        self.stdout.write(f"  Left alone: {counts['skip'] + counts['error']}")
        if reasons:
            self.stdout.write("  Reasons for leaving a row alone:")
            for reason, n in reasons.most_common():
                self.stdout.write(f"    {n:5d}  {reason}")
        if not apply_changes and counts["repair"]:
            self.stdout.write("")
            self.stdout.write("Dry run. Re-run with --apply to write these changes.")


def _qualifier(obj, name):
    try:
        value = obj.get_qualifier(name)
    except Exception:
        return None
    return str(value) if value is not None else None
