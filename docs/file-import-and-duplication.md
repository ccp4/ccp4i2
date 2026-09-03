# File Import and Duplication

How a file gets into a project, what that currently costs, and the one change
that would make the cost go away.

**Status:** the defects below are fixed (PR #347). The direction of travel —
a source store, so that one uploaded file can yield several representations —
is **designed but not started**, deliberately.

Companion documents:

- [PROJECT_STORE_AND_DB_LOCATION.md](PROJECT_STORE_AND_DB_LOCATION.md) — where
  the project directory and the database live, and how they are found.
- [PROJECT_RECOVERY.md](PROJECT_RECOVERY.md) — `DATABASE.db.xml`, the per-project
  snapshot this strand must not break.
- [MOVING_PROJECTS.md](MOVING_PROJECTS.md) — why absolute paths inside a project
  are load-bearing, which is why nothing here reaches across projects.

## The finding in one paragraph

Uploading a file copies it into the project with no notion of ever having seen
it before, and the record it writes cannot answer the question either:
`FileImport.checksum` is the checksum of the file *as it ended up in the
project*, so for anything derived — an MTZ split down to the columns one
parameter wanted — it is the checksum of the derivative, not of what the user
uploaded. Meanwhile the uploaded bytes are staged in `CCP4_DOWNLOADED_FILES`,
which nothing ever collects. The result is that fetching one reflection file for
F/SIGF and again for the free-R set — which is not a mistake, but the only way
to get both — stores the source twice, the conversion intermediate twice, and
two genuinely-different derived files. Nothing in the system can tell that those
four files came from one download.

## The two directories

The distinction already exists on disk and is worth stating, because the export
code already relies on it:

| | `CCP4_IMPORTED_FILES` | `CCP4_DOWNLOADED_FILES` |
|---|---|---|
| holds | representations | sources and conversion intermediates |
| has a DB row | a `File` (the only directory `File.path` can address) | nothing |
| offered in the file picker | yes | no |
| exported | per `File` row | wholesale, as a standard directory |
| ever collected | on replace | **never** |

That export asymmetry is not an accident to be tidied away — it is exactly the
split a source store wants, and it already round-trips.

Note that `CDataFile`'s path logic treats `CCP4_JOBS` and `CCP4_IMPORTED_FILES`
as *the* markers for "inside the project". A file under
`CCP4_DOWNLOADED_FILES` is therefore classified as **external** despite sitting
inside the project directory. Harmless while nothing points at those files; not
harmless once something does.

## What was fixed (PR #347)

1. **Files deleted that another job still needed.** Replacing the file a
   parameter holds unlinks the old one, and the guard saw only the container
   being edited — not a second *job* pointing at the same file, which
   "Browse the project hierarchy" produces routinely.
   `_is_referenced_outside_this_job` now answers from `FileUse` rows plus the
   params XML of jobs that have not run yet. When something else refers to the
   file, the bytes *and* the row are kept.
2. **Cross-project reuse.** `find_imported_file_by_checksum` searched every
   project; `File.path` resolves through `file.job.project.directory`, so a
   match elsewhere left a job pointing at bytes under a project it does not own.
   Now project-scoped, and it must stay that way.
3. **Failed uploads leaked their bytes.** Nothing owned a staged file until a
   `File` row existed, so a failure left it with nothing referencing it and
   nothing that would collect it. Two failed uploads of one 5.9 MB unmerged file
   left 11.8 MB.
4. **Unsplit uploads were stored twice.** An unmerged or selector-less MTZ
   derives nothing, so the staged and imported files were byte-identical twins.
   Moved rather than copied.
5. **`source_checksum`.** Taken from the bytes as uploaded, before any splitting
   or CIF conversion. Recorded but not acted on: a hash not taken at upload time
   cannot be recovered later, and this is the field every later step needs.

The endpoint now returns `duplicate_of` — earlier imports of the same source
bytes in this project, each flagged `interchangeable` — and acts on none of it.

**Same source does not mean interchangeable.** One reflection file imported for
F/SIGF and again for the free-R set yields two files that are both
`application/CCP4-mtz` but differ in content flag, and the picker would not
offer one for the other's parameter. Advising reuse there would be wrong, and
the fastest way to teach people to ignore the message. Only interchangeable
matches are surfaced. Matching is exact on type and content flag, which
under-advises (the picker also offers convertible content) — the safe direction
for a message whose only value is that it is trustworthy.

## What was deliberately not done

**Automatic reuse of a matching file.** It means two jobs sharing one `File`
row, and neither the delete path nor the export path is ready for that. The
unlink guard would have to become genuinely refcounted rather than
reference-*checked*, and `<importfile>` is one element per `File`.

**Chasing conversion intermediates on failure.** The cleanup removes only the
file `download_file` wrote. It cannot name the MTZ that `gemmi_convert_to_mtz`
writes beside a staged CIF, and guessing at neighbours by filename prefix risks
deleting a *concurrent* upload's staging — the same class of bug as item 1.

**Any garbage collection of `CCP4_DOWNLOADED_FILES`.** Sources are precisely
what the direction of travel wants to keep. Collecting them now would have to be
undone.

## Direction of travel: a source store

The model, in one sentence: **a `File` becomes a representation of a
`FileImport`**, and the byte-identical case is just the trivial representation.

On upload, hash the incoming bytes and compare against the `FileImport` rows of
this project. On a match:

- if a representation matching this parameter already exists, offer it;
- otherwise derive a new representation **from the retained source**, and attach
  it to the existing `FileImport`;
- either way, discard the incoming copy.

This is the only thing that fixes the F/SIGF-then-free-R case, because there the
second fetch is *required* under the current design — the user has done nothing
redundant, and no advisory can help them.

### What it commits to

- **`FileImport` stops being a satellite of one `File`.** Today it is
  `OneToOneField(File, primary_key=True)`. It must gain its own identity, with
  `File` pointing at it.
- **The source must be retained and addressable.** `FileImport` has no path
  today; `name` is only the client-side filename. Sources survive now purely
  because nothing deletes `CCP4_DOWNLOADED_FILES`. `FileImport` needs
  `rel_path` + `name`, allowed to point into *either* directory — for
  byte-identical imports the source *is* the representation, and storing a
  second copy to keep the model tidy would double disk on the common case.
- **Deletion becomes refcounted.** Unlinking a `File` must check whether a
  `FileImport` still needs those bytes for a future derivation.

### Blast radius

`FileImport` is read or written in seven non-test modules: `FileImportViewSet`,
`serializers.py`, `JobViewSet`, `async_db_handler`, `export_project`,
`import_i2xml`, `import_sqlite`.

Two are format-constrained and are the reason this has not been done casually:
the project export XML (`<importfile fileid=… sourcefilename=… checksum=…>`) and
the legacy `ImportFiles` table are **both keyed on `fileid`** — the Qt-era schema
was 1:1 too. The same tree is written per project as `DATABASE.db.xml`, which
exists for database recovery and for moving projects between machines. A reshape
can still round-trip — emit one element per `File`, duplicating the source
attributes, and re-converge on import by checksum — but that is a decision to
take deliberately, not to discover.

PR #347 keeps the one-to-one untouched for exactly this reason. Its only schema
change is one additive column, and its only serialization change is one
attribute that older readers ignore and whose absence reads as blank.

## On the measurements

Duplication was measured across the projects on one developer's machine. The
headline moved **21.3% → 27.5% → 11.9%** as the corpus was corrected — first for
omitting `CCP4_DOWNLOADED_FILES` entirely, then for missing a whole projects
root and for projects relocated by the `~/.ccp4i2x` home consolidation.

**Do not quote a percentage from this.** The corpus mixes scratch projects,
i2run test output and real work, and the figure swings by a factor of two
depending on which is included. Twelve of thirty-five projects had any
duplication at all; the median had none.

What *is* solid, because it was reproduced deliberately rather than inferred: a
new project driven through four uploads of one source plus two failures held
**14 MB in `CCP4_DOWNLOADED_FILES` against 1 MB in `CCP4_IMPORTED_FILES`**. And
in real projects, staging runs at 2–5x the imported data for several of them.

Before committing to the reshape, measure again on projects belonging to actual
crystallographers. The `duplicate_of` advisory now logs every hit server-side,
so that evidence accumulates without anyone collecting it by hand.

## Open items

- [ ] Measure duplication on real users' projects, not a developer's machine.
- [ ] Decide on the `FileImport` reshape; if yes, sequence it *after* an alpha,
      because it is a migration against installed user data.
- [ ] Garbage-collect `CCP4_DOWNLOADED_FILES` — but only once the source-store
      question is settled, since it decides what is garbage.
- [ ] Conversion intermediates: `gemmi_convert_to_mtz` writes an MTZ beside every
      staged CIF and nothing tracks it, on success or failure.
- [ ] `CDataFile` treats `CCP4_DOWNLOADED_FILES` as outside the project.
