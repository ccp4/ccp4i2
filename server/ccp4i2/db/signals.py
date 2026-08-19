"""Decide when a project's on-disk snapshot needs rewriting.

The rule is: anything the *user* authored that exists only in the database.
Everything about a project's files can be recovered from the project directory
itself; nothing about the user's judgement can. Job titles and comments,
evaluations of good and rejected, file annotations, tags and the project
description all vanish with the database unless they have been written out.

Driving this from signals rather than from calls at each endpoint means every
write path is covered -- REST, i2run, management commands, the job runner --
including paths that do not exist yet. DRF's ``partial_update`` saves without
``update_fields``, so which fields actually changed is worked out by comparing
against the stored row in ``pre_save``.

See :mod:`ccp4i2.db.project_snapshot` for how the writes are coalesced.
"""

import logging

from django.db.models.signals import m2m_changed, post_delete, post_save, pre_save
from django.dispatch import receiver

from .models import File, Job, Project, ProjectTag
from .project_snapshot import schedule_snapshot, update_registry

logger = logging.getLogger(f"ccp4i2:{__name__}")


#: Fields whose change means the snapshot is out of date. Everything here is
#: either user-authored or the job lifecycle, and none of it is reconstructible
#: from the project directory.
WATCHED_FIELDS = {
    Job: ("status", "title", "comments", "evaluation"),
    File: ("annotation", "sub_type", "content"),
    Project: ("name", "description", "directory", "last_job_number"),
}

#: Attribute used to carry the result of the pre_save comparison through to
#: post_save on the same instance.
_CHANGED = "_snapshot_changed_fields"


def _is_top_level(job: Job) -> bool:
    """Subjob transitions are frequent and their parent's snapshot covers them.

    A pipeline can run dozens of subjobs; snapshotting each time one starts or
    finishes would be constant churn for no gain, because the top-level job
    reaching a terminal state writes a snapshot that includes them all.

    ``parent`` is the truth here -- both create paths set it -- with the dotted
    job number as a fallback for rows that predate it or arrived by import.
    """
    if job.parent_id is not None:
        return False
    return "." not in str(job.number)


@receiver(pre_save, sender=Job)
@receiver(pre_save, sender=File)
@receiver(pre_save, sender=Project)
def record_changed_fields(sender, instance, **kwargs):
    """Note which watched fields differ from the stored row.

    Only for updates: on creation there is nothing to compare against, and the
    extra query would be paid on every row of a bulk import.
    """
    setattr(instance, _CHANGED, ())
    if instance.pk is None:
        return
    watched = WATCHED_FIELDS.get(sender)
    if not watched:
        return
    stored = sender.objects.filter(pk=instance.pk).only(*watched).first()
    if stored is None:
        return
    setattr(
        instance,
        _CHANGED,
        tuple(
            field
            for field in watched
            if getattr(stored, field, None) != getattr(instance, field, None)
        ),
    )


@receiver(post_save, sender=Job)
def job_saved(sender, instance, created, **kwargs):
    if created:
        # A new top-level job is a new piece of project history. Subjobs arrive
        # in bursts as a pipeline builds them and are covered by their parent.
        if _is_top_level(instance):
            schedule_snapshot(instance.project)
        return
    changed = getattr(instance, _CHANGED, ())
    if not changed:
        return
    # status is the lifecycle; the rest is the user's own writing, and worth
    # capturing on a subjob as well as a top-level one.
    if changed == ("status",) and not _is_top_level(instance):
        return
    schedule_snapshot(instance.project)


@receiver(post_save, sender=File)
def file_saved(sender, instance, created, **kwargs):
    # Files are created in bulk by the gleaner as a job finishes; that job's
    # own transition to a terminal state writes the snapshot that covers them.
    # What matters here is a later edit -- someone annotating a result.
    if created or not getattr(instance, _CHANGED, ()):
        return
    schedule_snapshot(instance.job.project if instance.job else None)


@receiver(post_save, sender=Project)
def project_saved(sender, instance, created, **kwargs):
    if created or getattr(instance, _CHANGED, ()):
        schedule_snapshot(instance)


@receiver(post_delete, sender=Job)
@receiver(post_delete, sender=File)
def row_deleted(sender, instance, **kwargs):
    # Deleting a project cascades to every job and file it owns, so by the time
    # these fire the project row may already be gone. Reaching for it then
    # raises; the project's own post_delete handles that case.
    try:
        project = (
            instance.project
            if isinstance(instance, Job)
            else (instance.job.project if instance.job_id else None)
        )
    except (Job.DoesNotExist, Project.DoesNotExist):
        return
    schedule_snapshot(project)


@receiver(post_delete, sender=Project)
def project_deleted(sender, instance, **kwargs):
    # The project's own snapshot goes with its directory; what needs correcting
    # is the registry that lists where projects are.
    update_registry()


@receiver(m2m_changed, sender=ProjectTag.projects.through)
def project_tags_changed(sender, instance, action, pk_set, **kwargs):
    if action not in ("post_add", "post_remove", "post_clear"):
        return
    if isinstance(instance, Project):
        schedule_snapshot(instance)
        return
    # instance is the tag; pk_set holds the projects it was applied to.
    for project in Project.objects.filter(pk__in=pk_set or []):
        schedule_snapshot(project)


@receiver(post_save, sender=ProjectTag)
def project_tag_saved(sender, instance, **kwargs):
    for project in instance.projects.all():
        schedule_snapshot(project)
