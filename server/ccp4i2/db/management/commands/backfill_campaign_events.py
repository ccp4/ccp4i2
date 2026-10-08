"""Rebuild the CampaignEvent projection from pandda_events receipts.

``CampaignEvent`` is a cache of what each finished receipt's ``params.xml``
records (docs/pandda-campaign-design.md, section 10.2). New receipts are
recorded as they finish; this fills in receipts that finished before the
table existed, and is the rebuild path after a project restore or import.

Idempotent: each receipt's rows are replaced, never added to.

Usage:
    python manage.py backfill_campaign_events                  # every receipt
    python manage.py backfill_campaign_events --group <uuid>   # one campaign's members
"""
from django.core.management.base import BaseCommand, CommandError

from ccp4i2.db.models import ProjectGroup
from ccp4i2.lib.campaign_events import backfill


class Command(BaseCommand):
    help = "Rebuild CampaignEvent rows from finished pandda_events receipts"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            "--group", type=str,
            help="Only the member projects of this campaign (its uuid; a name "
                 "or numeric id is accepted too)")

    def handle(self, *args, **options):
        group = None
        wanted = options.get("group")
        if wanted:
            group = self._find_group(wanted)
        counts = backfill(group)
        scope = f"campaign '{group.name}'" if group else "all projects"
        self.stdout.write(
            f"[OK] {counts['receipts']} receipt(s) read in {scope}; "
            f"{counts['events']} event row(s) written"
            + (f"; {counts['failed']} receipt(s) could not be read"
               if counts["failed"] else "")
        )

    @staticmethod
    def _find_group(wanted: str) -> ProjectGroup:
        queries = []
        try:
            import uuid
            queries.append({"uuid": uuid.UUID(wanted)})
        except ValueError:
            pass
        if wanted.isdigit():
            queries.append({"pk": int(wanted)})
        queries.append({"name": wanted})
        for query in queries:
            group = ProjectGroup.objects.filter(**query).first()
            if group is not None:
                return group
        raise CommandError(f"No campaign matches '{wanted}'")
