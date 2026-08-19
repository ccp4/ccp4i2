# Moving a project on disk

Users move projects: onto a bigger disk, into a tidier folder, off a laptop and
onto a share. Sometimes the path changes without them doing anything at all —
a mounted drive gets renamed by an administrator and every absolute path in the
project points at somewhere that no longer exists.

Qt-i2 could move a project but left jobs broken afterwards. This document
records where the breakage comes from and how the Django app handles it.

## The data model is nearly location-independent

Very little in CCP4i2 actually records where a project lives:

* `Project.directory` — the one absolute path in the database. `Job.directory`
  and `File.path` are both derived from it.
* `CDataFile` serialises as `project` (UUID) + `relPath` + `baseName` +
  `dbFileId`. Nothing absolute:

  ```xml
  <CPdbDataFile>
    <project>b18b4fa2450311f0bd4a56d135b7511f</project>
    <baseName>XYZIN_PARTIAL-coordinates.pdb</baseName>
    <relPath>CCP4_JOBS/job_2</relPath>
    <dbFileId>c069d452451311f0bd4a56d135b7511f</dbFileId>
  </CPdbDataFile>
  ```

`CDataFile.getFullPath()` is also mildly self-healing: an absolute `baseName`
that no longer exists is discarded in favour of `relPath` when one is set.

So a move is, to a first approximation, `mv` plus one `UPDATE`.

## Where absolute paths still leak in

Base classes are clean. The leakage is all in plugin-level code and in files
that programs write for themselves. Measured across twelve real projects:

| File | Written by | Re-read? |
|---|---|---|
| `program.xml` | the CCP4 program itself (aimless writes `<MTZmergedfilename>`) | by report classes |
| `params.xml`, `input_params.xml` | base serialiser — but only plugin-set scalars are absolute | clone, re-run |
| `*.scene.xml` | ccp4mg | yes — breaks the 3D view |
| `script.py`, `state*.py`, `*.scm` | Coot wrappers | yes — breaks Coot relaunch |
| `job.ccp4db.xml`, `DATABASE.db.xml` | legacy Qt-i2 | on import |
| `report.html`, `report_xml.xml` | report generator | cached, regenerable |
| `log.txt`, `stdout.txt`, `com.txt`, `*.log` | provenance | no |

### Plugin-level emitters, fixed at source

Three places were writing absolute paths into `params.xml`, where they then
propagated into every downstream job that inherited the parameter:

* `phaser_pipeline` set `KILLFILEPATH` to `<pipeline work dir>/INTERRUPT`. It no
  longer sets it at all: each phaser wrapper already defaults to
  `<its own work dir>/INTERRUPT`, which is also what
  `CPluginScript.isInterrupted()` watches.
* `aimless` stored `HKLOUT_BASENAME` as a full path. It now stores a bare
  basename and joins it with the work directory at the point of use. An
  absolute value inherited from an older job is stripped back to its basename,
  which incidentally fixes those jobs too.
* `phaser_EP_AUTO` built the "- original hand" annotation from
  `str(outputData.MAPOUT[0])` — the file object, whose `__str__` is its full
  path — instead of `str(...annotation)`. Corrected for `MAPOUT` and
  `LLGMAPOUT`.

Nothing can be done at source about `program.xml`: the program writes it.

## What a move does

`server/ccp4i2/db/move_project.py`.

1. **Preflight.** Destination must be absolute, must not exist, must have a
   writable parent, must not be inside the project. No job may be queued,
   running or running remotely — but a `PENDING` job is a task that has been
   created and not yet run, so those do not block. No other project may claim
   the destination.
2. **Move the bytes.** `os.rename` when the destination is on the same
   filesystem — atomic, instant, and no second copy of what can be a
   multi-gigabyte project. Cross-device falls back to `copytree`, and the
   source is only removed once the database has been updated, so an interrupted
   move always leaves one complete tree.
3. **Rebase.** Rewrite the old root to the new one in every text file that is
   not provenance, not a regenerable cache and not binary. Binaries are
   identified by sniffing for a NUL byte, not by trusting the extension.
4. **Delete the caches.** `report_xml.xml`, `report.html`, `report_tmp.html` —
   they regenerate on demand.
5. **Commit.** `Project.directory` is updated last. Any failure after step 2
   rolls back: the rewrites are inverted and the directory renamed back (or,
   cross-device, the partial copy discarded).

Logs, stdout, stderr and command scripts are deliberately left alone. They are
the record of what actually ran, and rewriting them would falsify it.

A journal file (`.ccp4i2-move.json`) sits at the destination root while the
rebase is in flight, so a tree left half-rebased by a crash can be recognised.

## Repairing paths without moving

A project may point somewhere that is neither its old nor its new home: a mount
renamed underneath it, an import from another machine, an earlier move. Every
operation reports these as `stale_roots` — a map of root to reference count —
and `repair_project_paths` rewrites any one of them onto the current directory,
in place, moving nothing.

## When the database itself points nowhere

The three operations above all assume the project directory is where
`Project.directory` says it is. A renamed drive breaks that assumption: the
database is intact, the files are intact, but the recorded path names a place
that no longer exists. Moving cannot help — there is nothing at the source to
move — and repairing cannot help either, because it rebases *onto* the recorded
directory.

`relocate_project` handles it: point a project at where its directory actually
is now, update the record, and rebase the paths inside from the old root to the
new one. Nothing is moved. The destination is checked for one of `CCP4_JOBS`,
`CCP4_IMPORTED_FILES`, `CCP4_COOT` or `CCP4_PROJECT_FILES` first, so a mistyped
or mis-picked folder is refused rather than rewritten through.

A drive rename invalidates *every* project at once, and re-pointing nineteen of
them one at a time is no kind of recovery path. `rebase_projects_root` takes the
old and new locations of the store, maps each project beneath the old root onto
its counterpart under the new one, and re-points the ones that are actually
there. Projects that cannot be found are reported and left alone — one failure
does not abandon the rest — and projects that were never on that drive are
untouched. `update_projects_dir_preference` follows the move in
`preferences.json` as well, so the next new project is not created back in the
location that has gone.

## API

| Endpoint | Purpose |
|---|---|
| `POST /projects/{id}/move/` | `{directory, dry_run?}` — relocate and rebase |
| `POST /projects/{id}/repair_paths/` | `{old_directory, dry_run?}` — rebase in place |
| `GET /projects/{id}/stale_roots/` | roots the project still refers to |
| `POST /projects/{id}/relocate/` | `{directory, dry_run?}` — re-point a project whose directory has gone |
| `GET /projects/missing_directories/` | projects whose recorded directory is not on disk |
| `POST /projects/rebase_root/` | `{old_root, new_root, dry_run?, update_preference?}` — the drive-rename case |

`PATCH`ing `Project.directory` directly is rejected by the serializer: it would
leave the record pointing at a directory still sitting somewhere else.

The frontend surfaces moving and repairing through `move-project-dialog.tsx`,
reached from the move icon on each row and card of the projects table, and
reconnection through `reconnect-projects-dialog.tsx`, reached from a banner that
appears when any project directory is missing. Both run the dry run
automatically, so the user sees what will change before anything is touched.
Projects that cannot be found are flagged in the list and do not offer a move,
which could not work until they are reconnected.

Both dialogs are desktop-only: the directories live on the machine running the
server, and choosing where they are needs a native folder picker.

## Known gap: importing a project zip

`import_ccp4_project_zip` accepts a `relocate_path` and updates
`Project.directory`, but does not rebase the paths inside the extracted files —
the same breakage this document describes, arriving by a different route.
Wiring the rebase into import is complicated by job renumbering: an imported
`job_3` may become `job_7`, so a plain root substitution is not enough; the
substitution list would have to be built from `import_i2xml_result["job_map"]`
as well as the root.
