"""Convert site tags into site evaluations.

Before the evaluation model, the only way to record what a dataset showed at
a site was to tag the project with the site's name -- and that tag was applied
to assert a hit. Nobody tagged a dataset at a site to say it was empty, or to
mark it as worth looking at later: the tag *was* the positive finding.

So these become HIT. Converting them to UNCLEAR would be the lossy choice: it
would discard a judgement somebody actually made and leave a campaign overview
claiming nothing is known about datasets whose hits were recorded.

The tag is then REMOVED. Keeping both would leave two overlapping records of
the same fact with nothing maintaining the older one -- the verdict control
writes evaluations, so a surviving tag would drift out of step and there would
be no way to tell which was current. The reverse migration puts the tags back,
so this is reversible rather than destructive.

The evaluations are marked as inferred rather than recorded -- evaluator
"migrated", plus a note -- because ProjectTag carries no author and no
timestamp, and a hit a person typed should stay distinguishable from one this
migration deduced. One query finds them all again if that reading turns out to
be wrong.

Tags that match no site are left exactly as they are. Under the old model the
association WAS the tag's text, so a site renamed since the tag was applied had
already lost its link before this migration ran -- there is nothing to recover,
which is precisely the orphaning that CampaignSite's stable ids now prevent for
good.
"""

from django.db import migrations

MIGRATED_EVALUATOR = "migrated"
MIGRATED_NOTE = (
    "Migrated from a site tag. Tagging a dataset at a site asserted a hit, so "
    "that is the verdict recorded here; it was inferred by migration 0025 "
    "rather than entered by a person."
)


def tags_to_evaluations(apps, schema_editor):
    ProjectTag = apps.get_model("ccp4i2", "ProjectTag")
    CampaignSite = apps.get_model("ccp4i2", "CampaignSite")
    SiteEvaluation = apps.get_model("ccp4i2", "SiteEvaluation")

    created = 0
    already = 0
    detached = 0
    unmatched = set()
    emptied = []

    for tag in ProjectTag.objects.all().iterator():
        for project in list(tag.projects.all()):
            # A tag is a SITE tag only if its text names a site of a campaign
            # the tagged project actually belongs to. Matching on text alone
            # would sweep up ordinary organisational tags that happen to share
            # a name with a site.
            #
            # Membership type is deliberately not filtered: a campaign's parent
            # can be tagged at a site like any other project, and excluding it
            # would silently drop that judgement.
            sites = list(
                CampaignSite.objects.filter(
                    group__memberships__project=project, name=tag.text
                ).distinct()
            )

            if not sites:
                # Only worth reporting if the text names a site somewhere --
                # otherwise it is an ordinary organisational tag.
                if CampaignSite.objects.filter(name=tag.text).exists():
                    unmatched.add((project.name, tag.text))
                continue

            for site in sites:
                if SiteEvaluation.objects.filter(
                    project=project, site=site
                ).exists():
                    already += 1
                    continue
                SiteEvaluation.objects.create(
                    project=project,
                    site=site,
                    verdict="hit",
                    evaluator=MIGRATED_EVALUATOR,
                    note=MIGRATED_NOTE,
                )
                created += 1

            # The evaluation now carries this fact, so the tag should not.
            tag.projects.remove(project)
            detached += 1

        # A tag that existed only to mark sites has nothing left to say.
        # Deleting it keeps the tag forest from filling with empty site names;
        # a tag that still labels other projects is left alone.
        tag.refresh_from_db()
        if not tag.projects.exists() and CampaignSite.objects.filter(
            name=tag.text
        ).exists() and not tag.children.exists():
            emptied.append(tag.text)
            tag.delete()

    if created or already or detached:
        print(
            f"\n  0025: {created} evaluation(s) created from site tags, "
            f"{already} pair(s) already evaluated, "
            f"{detached} tag link(s) removed."
        )
    if emptied:
        print(f"  0025: removed {len(emptied)} now-empty site tag(s).")
    for project_name, tag_text in sorted(unmatched):
        print(
            f"  0025: tag {tag_text!r} on {project_name!r} matches no site of "
            "its campaigns; left as it is."
        )


def evaluations_back_to_tags(apps, schema_editor):
    """Put the tags back and drop only the evaluations this migration created.

    Deleting by the "migrated" marker rather than by site, so verdicts a person
    has recorded since are not lost by a rollback. Each such evaluation names
    the site it came from, and the tag's text was that site's name, so the
    original tagging is reconstructible.
    """
    ProjectTag = apps.get_model("ccp4i2", "ProjectTag")
    SiteEvaluation = apps.get_model("ccp4i2", "SiteEvaluation")

    restored = 0
    for evaluation in SiteEvaluation.objects.filter(
        evaluator=MIGRATED_EVALUATOR
    ).select_related("site", "project"):
        tag, _ = ProjectTag.objects.get_or_create(
            text=evaluation.site.name,
            parent=None,
            defaults={"path": evaluation.site.name},
        )
        tag.projects.add(evaluation.project)
        restored += 1

    SiteEvaluation.objects.filter(evaluator=MIGRATED_EVALUATOR).delete()
    if restored:
        print(f"\n  0025 (reverse): restored {restored} site tag link(s).")


class Migration(migrations.Migration):

    dependencies = [
        ("ccp4i2", "0024_campaign_sites_and_evaluations"),
    ]

    operations = [
        migrations.RunPython(tags_to_evaluations, evaluations_back_to_tags),
    ]
