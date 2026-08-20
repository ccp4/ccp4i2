# Rebuilding the database from what is on disk

Everything CCP4i2 knows about a project's *files* can be recovered from the
project directory: the job directories are there, `params.xml` records what each
job ran with, the outputs are where they were written.

Everything the *user* wrote cannot. Annotations, job titles, comments, the
evaluation of one run as best and another as rejected, tags, the project
description — none of it has an on-disk artefact. It lives only in
`db.sqlite3`, and if that is lost or corrupted it is gone, with no way to
reconstruct it from anything.

Qt-era CCP4i2 knew this and wrote a per-project `DATABASE.db.xml` as jobs
started and finished, plus a list of known project directories outside the
database. The Django app had neither. This is that, restored.

## What gets written, and where

| File | Where | Holds |
|---|---|---|
| `DATABASE.db.xml` | project directory | everything the database knows about that project |
| `project_directories.json` | `~/.ccp4i2-django/` | where every known project is |

The registry matters as much as the snapshots. A snapshot inside a project
directory is only useful if something can find the project directory in the
first place, and with the database gone nothing else knows where the projects
are.

`DATABASE.db.xml` is the format the project export already produces and the
project import already reads — the one project zips carry and legacy project
directories already contain. Reusing it means one recovery path rather than
two. It also keys on UUIDs rather than on primary keys or on the current shape
of the models, which is what a file written to be read back *after something
went wrong* needs: primary keys are meaningless once the database has been
rebuilt, and a Django fixture would refuse to load after a migration changed
the models it was written against.

## When it gets written

Driven by database signals rather than by calls at each endpoint, so every
write path is covered — REST, `i2run`, management commands, the job runner —
including ones that do not exist yet.

The rule is: anything user-authored that exists only in the database.

| Model | Watched |
|---|---|
| `Job` | `status`, `title`, `comments`, `evaluation` |
| `File` | `annotation`, `sub_type`, `content` |
| `Project` | `name`, `description`, `directory`, `last_job_number` |
| `ProjectTag` | tag text and which projects carry it |

Two deliberate exclusions:

- **Subjob status changes.** A pipeline runs dozens of subjobs; the top-level
  job reaching a terminal state writes a snapshot covering them all. A comment
  or evaluation on a subjob *is* captured — status churn is noise, what the
  user wrote is not.
- **File creation.** The gleaner drops many files as a job finishes, and that
  job's own transition covers them. A later edit — someone annotating a result
  — is what triggers a write.

DRF's `partial_update` saves without `update_fields`, so which fields actually
changed is worked out by comparing against the stored row in `pre_save`.

## Coalescing, and why there is no timer

One request can touch many rows. Snapshots are registered with
`transaction.on_commit` and de-duplicated per transaction, so a request that
changes thirty files still writes once.

That is coalescing **by transaction, not by wall-clock**. There is deliberately
no timer anywhere here, because there is nowhere for one to live: the desktop
app runs uvicorn with two workers, and `i2run` and the job runner are processes
that exit as soon as they are done. An in-process timer would be duplicated
across the workers and would never fire at all in the CLI.

The de-duplication set is Django's own `run_on_commit` list rather than one
kept alongside it, because Django clears that list on rollback as well as on
commit. A separate set would keep an entry for a transaction that rolled back,
and the next change to that project would be silently skipped — a safety net
with a hole in it that nothing would ever report.

Writes are atomic (temp file, `os.replace`), so two workers writing at once
cannot produce a partial file. Both read the same database, so last-writer-wins
is harmless.

Failing to write a snapshot never propagates. A job must not be reported as
failed because its directory was read-only: a safety net that breaks the thing
it protects is worse than no safety net.

## Retrofitting projects that predate this

Signals only fire when something changes, and a finished, dormant project is
both the kind most worth protecting and the kind nothing will ever touch again.
So there has to be a way to write them all out:

```bash
ccp4-python manage.py snapshot_projects --status   # who is protected
ccp4-python manage.py snapshot_projects            # write the missing ones
ccp4-python manage.py snapshot_projects --all      # rewrite everything
```

The desktop app runs `snapshot_projects` at launch, just after `migrate`. In
its default missing-only mode it is a no-op once every project is covered, so
it is cheap to run every time, and it is never fatal.

### Legacy snapshots are preserved, not overwritten

A project directory carried over from Qt-era CCP4i2 already has a
`DATABASE.db.xml`, and it may describe jobs or annotations that were never
imported into this database. The first time one is about to be overwritten, the
original is copied to `DATABASE.db.xml.pre-ccp4i2-django` — once, so our own
output can never creep over the real backup.

Snapshots we wrote are recognised by a `<generator>ccp4i2-django</generator>`
element in the header. Recognition matches the whole element, not the marker on
its own: the first real project this was tried on referenced a directory called
`ccp4i2-djangodbapi`, and a substring search mistook it for one of ours.

## Suspending

Bulk work — importing a project zip, migrating a legacy database — would
otherwise snapshot after every row:

```python
from ccp4i2.db import project_snapshot

with project_snapshot.suspended():
    ...  # hundreds of rows
project_snapshot.write_snapshot(project)
```

`import_ccp4_project_zip` already does this, snapshotting once at the end when
the files are actually on disk. Setting `CCP4I2_DISABLE_PROJECT_SNAPSHOT=1`
disables writing entirely.

## Rebuilding

Restoring is the other half, and it is deliberately *not* the same operation as
importing, because the two have opposite invariants.

**Import** brings a project in alongside an existing world. Another copy may
already be present, and job numbers may collide with jobs that already exist,
so renumbering is correct — `import_ccp4_project_zip` renames the job
directories as it extracts them so disk and database still agree afterwards.

**Restore** is the reverse. The directory is the truth and the database is what
is being repaired. Nothing may be renumbered, because the job directories are
already on disk under the names they have and nothing is going to move them. A
clash is therefore an error to report, not something to work around, and the
post-condition is checked rather than assumed: after a restore, every job the
database records must have a directory on disk. `verify_restored` checks
exactly that, which catches renumbering, a mis-rooted project and a snapshot
written for a different tree — none of which would surface as an error during
the import itself.

```bash
manage.py restore_projects --registry              # everything the registry knows
manage.py restore_projects --scan ~/CCP4I2_PROJECTS  # registry lost as well
manage.py restore_projects --directory <one project>
manage.py restore_projects --registry --dry-run    # look first
```

A project already in the database is left alone unless `--replace` is given:
overwriting a live project with an older snapshot would lose work rather than
recover it. One project failing does not abandon the rest.

If a project directory has moved since its snapshot was written, the directory
wins — the snapshot is re-rooted onto it, and the absolute paths inside the
project's files are rebased to match, reusing the machinery from
`docs/MOVING_PROJECTS.md`.

`GET /projects/restorable/` (optionally `?scan=PATH`) previews what could be
recovered without touching anything; `POST /projects/restore/` does it. The
desktop app exposes both through a Recover button on the projects toolbar.

## What this does not do

Inferring a project from its directory when there is **no** snapshot at all is
still missing, and is necessarily lossy. See the notes on
`utils/reconstructDBFromXML.py`, whose own TODO header is an honest account of
how far inference gets you — it recovered 553 of 577 files on the project it
was measured against, and *over-generated* file-use rows, 269 against a true
209.
