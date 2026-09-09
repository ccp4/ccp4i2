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

`captureImportProvenance` is a per-browser UI preference. It is **on by
default** (the Qt-era behaviour). Two ways to change it, neither needing a
settings pane:

- The prompt itself carries a **"Don't ask again"** checkbox. Ticking it and
  finishing (Save or Skip) turns the preference off from that first exposure
  on — this import's note is still taken, and no further imports prompt.
- The **View** menu toggles it either way ("Ask for Import Provenance" /
  "Don't Ask for Import Provenance"), so a user who opted out can opt back in.

The preference is read at each import (`readUiPreference`), so a change from
either route takes effect on the very next import.

## Coverage

Only four components call `uploadFileParam` directly; two of them own every
user-initiated **local-file** import, and both are wired:

| Uploader | Wired? | Covers |
|---|---|---|
| `csimpledatafile` (`CSimpleDataFileElement`) | ✅ | generic data files, **coordinates** (`cpdbdatafile` renders it), sequences, dictionaries, TLS, and the raw MTZ pick in `import_merged` (`CGenericReflDataFile` HKLIN) |
| `cminimtzdatafile` (`CMiniMtzDataFileElement`) | ✅ | **MTZ** obs/map/phases (`CObsDataFile` etc.), and **free-R** (`cfreerfile` wraps it); one prompt on the pick, applied to both the F/SIGF and the free-R-sibling upload |
| `fetch-file-for-param` (fetch from the internet / PDB) | — | *deferred*: the source is a URL/accession, already self-describing; a future touch could auto-record it as the note without prompting |

`import_merged` needs no wiring of its own: the user picks HKLIN through
`csimpledatafile` (which prompts and stores the note on that file), and the
subsequent split into `HKLIN_OBS` is a *derived* upload — see the rule below.
Task interfaces that render standard file elements — **`splitMtz`**, the
`Import*` family — are likewise covered through the element uploaders.

### One prompt per pick — no dedup

The prompt fires at the single moment the user picks a file from disk. A
monolithic MTZ split into F/SIGF and free-R, or `import_merged` re-uploading a
split of the file just picked, produces *further* `uploadFileParam` calls —
but those are **derived, not user picks**, so they never call
`requestImportProvenance` and never prompt. The note is captured once, at the
pick, and stored on the file the user chose. Because there is only ever one
call per pick, there is no cross-call dedup and no timing window to reason
about. (`cminimtzdatafile` captures the note once and passes it to both the
F/SIGF and free-R uploads explicitly.)

### Adding it to a new user-pick uploader

```ts
const { requestImportProvenance } = useImportProvenance();
// ... in the user-pick handler, before uploadFileParam:
const provenance = await requestImportProvenance(file.name);
await uploadFileParam({ /* ... */, description: provenance ?? undefined });
```

Call it **only** for a genuine from-disk pick, and pass the note to every
upload that pick produces (including derived splits). Never call it from a
programmatic upload that has no user behind it.

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
