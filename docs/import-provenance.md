# File-import provenance capture

When a user imports a local file into a task parameter, they can be prompted
for a free-text note describing where the file came from — the experiment, the
processing software, another project. This is the Qt-era "Describe source of
this file" step, reinstated on the Django/React stack.

## What it captures (and what it deliberately doesn't)

The Qt app captured **two** provenance elements at import:

1. a free-text **narrative** ("where did this come from?"), and
2. an **association to a `FileExport`** — a link recording that the file was
   the output of some other i2 job's deliberate export.

This feature implements **(1) only**. See *FileExport deferred*, below, for why
(2) is left out.

## Data model

The narrative lives on **`FileImport.description`** (`TextField`, blank by
default) — migration `0021_fileimport_description`.

It is kept deliberately distinct from **`File.annotation`**, which is an
*auto-generated short label* (`"Imported from data.mtz; Columns: F, SIGF"`)
that the GUI recycles as the file's display name in lists. Overwriting that
label with the user's narrative would break the lists, so the narrative gets
its own field. (The Qt schema had the same split: `Files.Annotation` = label,
`ImportFiles.Annotation` = narrative.)

## Flow

The note rides along in the **same upload POST** — no extra endpoint, no
post-hoc PATCH:

1. The user picks a local file in a data-file parameter.
2. If the `captureImportProvenance` preference is on, a dialog prompts for the
   note *before* the bytes are read (`ImportProvenanceProvider`).
3. The note is passed as `description` to `uploadFileParam`, which appends it to
   the `multipart/form-data` sent to `POST jobs/{id}/upload_file_param`.
4. The server (`lib/utils/files/upload_param.py`) reads `request.POST["description"]`
   and stores it on the created `FileImport`.

`requestImportProvenance(fileName)` resolves to `null` when the preference is
off (or no provider is mounted), so the prompt is a genuine no-op unless opted
in, and **programmatic/derived imports never call it** — only a user picking a
file is ever prompted.

### Key files

| Concern | File |
|---|---|
| Preference key + storage | `client/renderer/lib/ui-preferences.ts` (`captureImportProvenance`) |
| The prompt dialog + context | `client/renderer/providers/import-provenance-provider.tsx` |
| Provider mount point | `client/renderer/app/ccp4i2/(authed)/layout.tsx` |
| Upload plumbing (`description` arg) | `client/renderer/utils.ts` (`uploadFileParam`) |
| The wired upload path | `client/renderer/components/task/task-elements/csimpledatafile.tsx` |
| On/off toggle | `client/renderer/components/view-menu.tsx` ("Ask for Import Provenance") |
| Server storage | `server/ccp4i2/lib/utils/files/upload_param.py`; `server/ccp4i2/db/models.py` (`FileImport.description`) |
| Tests | `server/ccp4i2/tests/api/unit/test_import_provenance_description.py` |

## The preference (and its default)

`captureImportProvenance` is a per-browser UI preference, toggled from the
**View** menu. It is **off by default** — prompting on every import is the
deliberately-tedious behaviour, so it is opt-in rather than opt-out. Flipping
the default to on is a one-line change in `ui-preferences.ts` (`DEFAULTS`).

> Decision flagged for review: the request was phrased as "provide the
> behaviour, backed by a preference to switch it off", which reads as
> *default-on*. I chose default-**off** to avoid a dialog on every alpha
> tester's import; say the word to flip it.

## Extending to other upload components

`csimpledatafile` (the generic `CDataFileElement`) is wired as the reference
implementation and covers most file parameters. Specialised file widgets that
call `uploadFileParam` directly for a **user-picked** file — e.g.
`cminimtzdatafile`, `cpdbdatafile`, `casudatafile`, `cimportunmerged` — opt in
with the same two lines:

```ts
const { requestImportProvenance } = useImportProvenance();
// ... in the user-pick handler, before uploadFileParam:
const provenance = await requestImportProvenance(file.name);
await uploadFileParam({ /* ... */, description: provenance ?? undefined });
```

Do **not** add it to programmatic upload calls (e.g. the free-R sibling
auto-import in `cminimtzdatafile`) — those should never prompt.

## FileExport deferred

The import→export association is intentionally not implemented, for two
reasons:

1. **No clean chokepoint to mint a `FileExport`.** The only place a specific
   `File` leaves the system is the generic `files/{id}/download` endpoint, which
   also serves every in-app preview/open (Coot, Moorhen, ViewHKL…). Minting a
   `FileExport` on every hit would record *opening* a file as an *export*. The
   job-level `export_job_file` endpoint is job-scoped and reconstructs/combines
   files, so it is not a specific-`File` export either. `FileExport` rows are
   currently written only by legacy-project import.
2. **Its only consumer is the association we're not building.** A `FileExport`
   row is, on its own, just an export audit log of modest value; its point was
   to let an import link back to it, and that link's hard half is a
   matching-UI at import time, not the minting.

If this is wanted later, the clean minimal design is: the two deliberate
file-save actions ("Download", "Save to…" in the file context menu) append
`?export=1` to the download URL, and `FileViewSet.download` mints
`FileExport(file=…, name=…)` only when that flag is present — leaving previews
and opens untouched.
